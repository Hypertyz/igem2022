# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 00:20:56 2022

@author: HoangPham
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# writing to file
# Using readlines()
file1 = open('orf_coding_all_R64-3-1_20210421.fasta', 'r')
Lines = file1.readlines()


out = []
# Strips the newline character
for idx, line in enumerate(Lines):
    if ">" in line:
        out.append(line[1:])

f = open("gene_info.csv", "a")
f.writelines(out)
f.close()
