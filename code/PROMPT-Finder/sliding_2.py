# -*- codeing = utf-8 -*-
# @Time :  1:53 PM
# @Author : Berial
# @File: sliding_2.py
# @Software: PyCharm

# this version is add fdr to make result more accurate

##work for antisense chromatin RNA-seq
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.multitest import fdrcorrection as fdr
import matplotlib.pyplot as plt
import argparse
import concurrent.futures
import random
import time
import os

parser = argparse.ArgumentParser(description='process sliding window for chromatin RNA-seq')
parser.add_argument("-f", "--forward", type=str, metavar="forward", required=True, help="chromatin forward strand bedgrph file")
parser.add_argument("-r", "--reverse", type=str, metavar="reverse", required=True, help="chromatin reverse strand bedgrph file")
parser.add_argument("-b", "--binsize", type=int, metavar="binsize", default=200, help="binsize for sliding window default = 200bp")
parser.add_argument("-s", "--sliding", type=int, metavar="sliding", default=10, help="sliding length for sliding window. default = 10bp")
parser.add_argument("-o", "--output", type=str, metavar="output", default="result.bed", help="the output file name. default = result.bed")
parser.add_argument("-c", "--cutoff", type=float, metavar="cutoff", default=0.05, help="the pvalue cutoff of probability distribution. default = 0.05")
parser.add_argument("-p", "--threads", type=int, metavar="threads", default=23, help="the maximum threads allow to use. default = 23")
parser.add_argument("-sp", "--species", type=str, metavar="species", default="hg38", help="the specie of data")

args = parser.parse_args()
_forward_file = args.forward
_reverse_file = args.reverse
_bin = args.binsize
_sliding = args.sliding
_output_file = args.output.strip(".bed")
_threads = args.threads

start = time.perf_counter()
##read genome size
with open("/reference/{}/sizes.genome".format(args.species), "r") as f1:
    _size = {}
    for line in f1:
        chr, len = line.strip().split("\t")
        if "_" not in chr:
            if int(len) > 2000000:
                _size[chr] = int(len)


def add_key(dict, k):
    if k not in dict.keys():
        dict[k] = ""


# generate random files
def rando_bin(_chr, _strand):
    print("start {} {}".format(_chr, _strand))
    # read known site
    with open("/reference/{}/RefSeq_all.txt".format(args.species), "r") as f2:
        _rna = {}
        for line in f2:
            if "{}\t".format(_chr) in line:
                sep = line.strip().split("\t")
                chr, strand, start, end = sep[2:6]
                end = int(end)
                start = int(start)
                if _strand == "+":
                    end = end + 10000
                    ASstart = start - 20000
                    ASend = start
                elif _strand == "-":
                    start = start - 10000
                    ASstart = end
                    ASend = end + 20000
                if strand == _strand:
                    for i in range(start, end):
                        add_key(_rna, i)
                else:
                    for i in range(ASstart, ASend):
                        add_key(_rna, i)
    _rna = _rna.keys()
    ## read bg file
    _bg = {}
    if _strand == "+":
        file = _forward_file
    else:
        file = _reverse_file
    with open(file, "r") as f1:
        for line in f1:
            if line.startswith(_chr + "\t"):
                chr, start, end, count = line.strip().split("\t")
                start = int(start)
                end = int(end)
                count = int(count)
                if chr == _chr and count != 0:
                    for i in range(start, end):
                        _bg[i] = count
                if chr != _chr:
                    break
    # random results
    _random = []
    size = _size[_chr]
    chr = _chr
    # get random bins
    for _ in range(0, 10000):
        total = 0
        _bin_start = random.randint(0, size - _bin)
        _bin_end = _bin_start + _bin
        while _bin_start in _rna or _bin_end in _rna:
            _bin_start = random.randint(0, size - _bin)
            _bin_end = _bin_start + _bin
        for i in range(_bin_start, _bin_end):
            total += _bg.get(i, 0)
        _random.append(total)
    # calculate ecdf
    ecdf = ECDF(_random)
    # generate pvalue_list bed_list to record
    pvalue_list = []
    bed_list = []
    # generate pvalue list in genome wide
    chr = _chr
    start = 0
    end = _size[chr]
    _total = []
    for start1 in range(start, (end - _bin), _sliding):
        if not _total:
            for i in range(start1, (start1 + _bin)):
                _total.append(_bg.get(i, 0))
        else:
            del _total[:_sliding]
            for i in range((start1 + _bin - _sliding), (start1 + _bin)):
                _total.append(_bg.get(i, 0))
        if start1 + _bin >= end:
            print("finished {} {}".format(_chr, _strand))
            break
        if 0 in _total:
            continue
        total = sum(_total)
        pvalue = 1 - ecdf(total)
        bed = "{},{},{},{}".format(chr, start1, start1 + _bin, _strand)
        pvalue_list.append(pvalue)
        bed_list.append(bed)
    padj_list = fdr(pvalue_list)[1]
    with open(_output_file + "temp", "a") as fo:
        for padj, bed in zip(padj_list, bed_list):
            if padj > args.cutoff:
                continue
            chr, _start, _end, _strand = bed.split(",")
            fo.write("\t".join([chr, str(_start), str(_end), ".", ".", _strand]) + "\n")


def main(i):
    rando_bin(i, "+")
    rando_bin(i, "-")


if __name__ == "__main__":
    chr = []
    for i in range(1, 23):
        chr.append("chr" + str(i))
    chr.append("chrX")
    with concurrent.futures.ProcessPoolExecutor(max_workers=_threads) as executor:
        executor.map(main, chr)
    os.system("sort -k1,1 -k2,2n {}>{}".format(_output_file + "temp", _output_file + "temp1"))
    os.system("bedtools merge -s -c 6,6,6 -o distinct -i {} > {}".format(_output_file + "temp1", args.output))
    os.system("rm {}temp {}temp1".format(_output_file, _output_file))
    end = time.perf_counter()
    print("all done")
    time = end - start
    hour = int(time / 3600)
    min = int(time / 60) - hour * 60
    seconds = int(time - hour * 3600 - min * 60)
    print("time: {}hour(s) {}minute(s) {}second(s)".format(hour, min, seconds))
