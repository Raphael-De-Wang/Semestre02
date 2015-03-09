#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *

fastaName           = "AT_arc_metatrans.filtered.fasta.cleanup"
step2Name           = "AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step"
gCAIsNivExprFName   = "gCAIsNivExpr.txt"

fastaDict = readFASTA(fastaName)
step2List = read2step(step2Name)

with open(gCAIsNivExprFName) as fb:
    for line in fb:
        larray = line.split()
        ID     = larray[0]
        domain = ID.split('_')[-1]
        gCAIs  = float(larray[1])
        if gCAIs < 0.95:
            break
        writeFileName= fastaName + '.' + domain
        trierFastaByDomain(domain,fastaDict,step2List,writeFileName,Step2FormatSepArchs)

