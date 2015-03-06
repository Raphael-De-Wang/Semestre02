#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *

WORKPATH="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/"

def Step2FormatSepBestDomain(line):
    score = float(line[0])
    gId   = line[4]
    gName,ARC,RF = gId.split("__")
    domain = line[5]
    RF = RF.split("_")[-1]
    reverse,begin,end = indiceBeginEnd(RF,int(line[2]),int(line[3]))
    gId   = "%s_%s_%d_%d_%s"%(gName,RF,begin,end,domain)
    return (score,gName,gId,ARC,RF,reverse,begin,end,"")

def Step2FormatSepArchs(line):
    score = float(line[0])
    gId   = line[3]
    domain= line[4]
    gName,ARC,RF = gId.split("__")
    RF = RF.split("_")[-1]
    reverse,begin,end = indiceBeginEnd(RF,int(line[1]),int(line[2]))
    gId   = "%s_%s_%d_%d_%s"%(gName,RF,begin,end,domain)
    return (score,gName,gId,ARC,RF,reverse,begin,end,"")

fastaDict = readFASTA(WORKPATH+"Scripts/AT_arc_metatrans.filtered.fasta")
# step2M = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
# cleanUpFasta(WORKPATH+"Scripts/AT_arc_metatrans.filtered.fasta.cleanup",fastaDict,step2M,Step2FormatSepArchs)
step2M = read2step("best.domains.2step")
cleanUpFasta(WORKPATH+"Scripts/AT_arc_metatrans.filtered.fasta.cleanup",fastaDict,step2M,Step2FormatSepBestDomain,seuil=1e-8)
