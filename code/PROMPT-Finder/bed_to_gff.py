# -*- codeing = utf-8 -*-
# @Time : 2021/5/9 11:46 AM
# @Author : Berial
# @File: bed_to_gff.py
# @Software: PyCharm
##write bed to gff file
import sys
import os
if sys.argv[0] == "-h" or ".bed" not in sys.argv[1]:
    print("work to transform bedfile to gff only with trancsripts\n"
          "usage:\n"
          "python bed_to_bed_transcripts.py input.bed output.gff\n")
    exit(0)

os.system("sort -k1,1 -k2,2n {} > temp.bed".format(sys.argv[1]))
f1 = open("temp.bed", "r")
fo = open(sys.argv[2], "w")
fo.write("# creat from bed file and only have transcriptsID with out exon information\n")
for line in f1:
    chr, start, end, name, score, strand = line.strip().split("\t")[0:6]
    if int(end) > int(start):
        if "_" not in chr:
            fo.write("\t".join([chr, "Berial", "transcript", start, end, ".", strand, ".", "transcript_id="+name]) + "\n")
f1.close()
fo.close()
os.system("rm temp.bed")