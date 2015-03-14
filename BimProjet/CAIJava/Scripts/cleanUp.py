#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *

WORKPATH="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/"

fastaDict = readFASTA(WORKPATH+"Scripts/AT_arc_metatrans.filtered.fasta")
step2M = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
cleanUpFasta(WORKPATH+"Scripts/AT_arc_metatrans.filtered.fasta.cleanup",fastaDict,step2M,Step2FormatSepArchs)

# step2M = read2step("best.domains.2step")
# cleanUpFasta(WORKPATH+"Scripts/resulats/AT_arc_metatrans.filtered.fasta.cleanup",fastaDict,step2M,Step2FormatSepBestDomain,seuil=1e-8)
