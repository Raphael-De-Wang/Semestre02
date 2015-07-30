#!env python
import sys
import numpy as np

# tgt = 'outtxt'
tgt = sys.argv[1]

vec   = {}
count = 0
nflag = True
flag  = False
flag2 = False
handle = open(tgt)

'''
for line in handle:
    if nflag is True:
        nflag = False
        vname = line.strip().split('/')[-1].split('.')[0]
    if flag and count < 8:
        if '(1%)' in line:
            continue
        sp = line.split()
        for s in sp:
            key,value = s.split(':')
            vec[key] = value
        count += 1
    if 'Iteration 14' in line:
        flag = True
'''

for line in handle:
    if nflag is True:
        nflag = False
        vname = line.strip().split('/')[-1].split('.')[0]
    if flag and flag2 and count < 8:
        sp = line.split()
        for s in sp:
            key,value = s.split(':')
            vec[key] = value
        count += 1
    if 'CodonUsageIt14' in line:
        flag = True
    if flag and len(line.strip()) == 0 :
        flag2 = True

handle.close()

keys = ['tat', 'tgt', 'ggt', 'tct', 'ttt', 'tgc', 'tag', 'taa', 'tac', 'ttc', 'tcg', 'tta', 'ttg', 'tcc', 'tca', 'gca', 'gta', 'gcc', 'gtc', 'gcg', 'gtg', 'cgt', 'gtt', 'gct', 'gat', 'ctt', 'cct', 'cga', 'cgc', 'ctc', 'aca', 'cgg', 'ggg', 'gga', 'ggc', 'gag', 'acg', 'gac', 'ccg', 'gaa', 'acc', 'atg', 'aag', 'aaa', 'atc', 'aac', 'ata', 'agg', 'cag', 'agc', 'aga', 'cat', 'aat', 'att', 'ctg', 'cta', 'act', 'cac', 'tga', 'caa', 'agt', 'cca', 'ccc', 'tgg']

print '"'+vname+'" '+' '.join([ vec[k].replace(',','.') for k in keys ])

