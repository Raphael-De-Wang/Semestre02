#!env python

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

import numpy as np
import matplotlib.pyplot as plt

count = 0.

aci_stat_dict = {"tat":0, "tgt":0, "ggt":0, "tct":0, "ttt":0, "tgc":0, "tag":0, "taa":0,
                 "tac":0, "ttc":0, "tcg":0, "tta":0, "ttg":0, "tcc":0, "tca":0, "gca":0,
                 "gta":0, "gcc":0, "gtc":0, "gcg":0, "gtg":0, "cgt":0, "gtt":0, "gct":0,
                 "gat":0, "ctt":0, "cct":0, "cga":0, "cgc":0, "ctc":0, "aca":0, "cgg":0,
                 "ggg":0, "gga":0, "ggc":0, "gag":0, "acg":0, "gac":0, "ccg":0, "gaa":0,
                 "acc":0, "atg":0, "aag":0, "aaa":0, "atc":0, "aac":0, "ata":0, "agg":0,
                 "cag":0, "agc":0, "aga":0, "cat":0, "aat":0, "att":0, "ctg":0, "cta":0,
                 "act":0, "cac":0, "tga":0, "caa":0, "agt":0, "cca":0, "ccc":0, "tgg":0, "nnn":0 }

def acid_amino_statis (seq, aci_stat_dict):
    for fr in range(0,len(seq),3):
        try:
            aci_stat_dict[seq[fr:fr+3]] += 1
        except:
            aci_stat_dict["nnn"] += 1
    return len(seq)/3

def acid_amino_statis_job (seq, aci_stat_dict):
    with open(sys.argv[1]) as ihandle:
        for record in SeqIO.parse(ihandle, 'fasta'):
            count += acid_amino_statis (str(record.seq.lower()), aci_stat_dict)
    return count

acid_amino_list = ["tat", "tgt", "ggt", "tct", "ttt", "tgc", "tag", "taa", "tac", "ttc", "tcg", "tta", "ttg", "tcc", "tca", "gca", "gta", "gcc", "gtc", "gcg", "gtg", "cgt", "gtt", "gct", "gat", "ctt", "cct", "cga", "cgc", "ctc", "aca", "cgg", "ggg", "gga", "ggc", "gag", "acg", "gac", "ccg", "gaa", "acc", "atg", "aag", "aaa", "atc", "aac", "ata", "agg", "cag", "agc", "aga", "cat", "aat", "att", "ctg", "cta", "act", "cac", "tga", "caa", "agt", "cca", "ccc", "tgg" ]

sign_list = [ float(s) for s in "0.023 0.096 0.616 0.156 0.023 1.000 1.000 0.235 1.000 1.000 0.034 0.005 0.110 1.000 0.018 0.015 0.003 1.000 1.000 0.019 0.345 0.508 0.207 0.270 0.258 0.210 0.122 0.013 0.548 1.000 0.010 0.016 0.043 0.525 1.000 1.000 0.022 1.000 0.066 0.136 1.000 1.000 1.000 0.029 1.000 1.000 0.000 0.799 1.000 0.143 1.000 0.099 0.013 0.102 0.408 0.003 0.114 1.000 0.000 0.075 0.016 0.076 1.000 1.000".split(" ") ]

def print_acid_amino_percent(acid_amino_list, aci_stat_dict, count) :
    outstr = ""
    for aa in acid_amino_list :
        outstr += "%f "%(aci_stat_dict[aa]/count)
    print outstr
    
def plot_acid_amino_percent_contre_signature(acid_amino_list, sign_list, aci_stat_list) :
    fig = plt.figure()
    plt.xticks(range(64),acid_amino_list)
    for ind,sa in enumerate(zip(sign_list, aci_stat_list)):
        plt.bar(ind, sa[0], color='r', width=0.5)
        plt.bar(ind, 20.*sa[1]/count, color='b', width=0.3)
    plt.grid(True)
    plt.show()

def correcoef_sign_acid_count():    
    # aci_stat_list = [ aci_stat_dict[aa] for aa in acid_amino_list ]
    # plot_acid_amino_percent_contre_signature(acid_amino_list, sign_list, aci_stat_list)
    # coefficient
    print np.corrcoef(sign_list, aci_stat_list)

def correcoef_cais_1_2(inum):
    orf1 = []
    orf2 = []
    with open("/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/domain_cais_1_2_len%d.lst"%inum) as handle :
        for line in handle:
            (lable,v1,v2) = line.split()
            orf1.append(float(v1))
            orf2.append(float(v2))
    print np.corrcoef(orf1,orf2)

for i in range(1,9):
    correcoef_cais_1_2(i)
