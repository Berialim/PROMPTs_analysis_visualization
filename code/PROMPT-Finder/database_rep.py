# -*- codeing = utf-8 -*-
# @Time : 2021/8/19 7:46 PM
# @Author : Berial
# @File: database_rep.py
# @Software: PyCharm
# add a col for rep1 rep2
import os
import shutil
sample = {}
repeat = {}
with open("database.csv") as f:
    for line in f:
        if line.startswith("name"):
            continue
        name, condition = line.strip().split(",")
        if condition in sample.keys():
            sample[condition] += [[name, "rep_" + str(len(sample[condition]) + 1)]]
        else:
            sample[condition] = [[name, "rep_1"]]
for k, v in sample.items():
    for i in v:
        name, rep = i
        if rep in repeat.keys():
            repeat[rep] += [[name, k]]
        else:
            repeat[rep] = [[name, k]]
for k, v in repeat.items():
    for i in v:
        name, condition = i
        if not os.path.exists(k):
            os.mkdir(k)
        shutil.copyfile(name + ".bam", k + "/" + condition + ".bam")
        # shutil.move(name + ".bam", k)
        # shutil.move(name + ".bam.bai", k)
        # os.chdir(k)
        # os.renames(name + ".bam", condition + ".bam")
        # os.renames(name + ".bam.bai", condition + ".bam.bai")
        # os.chdir("..")