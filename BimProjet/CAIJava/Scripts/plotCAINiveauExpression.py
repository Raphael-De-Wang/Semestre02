#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *

# loading
step2M   = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
print step2M[binSearch(step2M,step2GnameComp,"GG7SD3401DYJ3N")]
print nivExprDomain(step2M,4)
exit()
CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
CAI_dict = readgCAIs()

CLS_sorted = np.array(sortDictByValue(CLS_dict))
# gCAIs_sorted = np.array(sortDictByValue(CAI_dict))


'''
# strategy 1 mean ou max of duplications
gName   = []
nivExpr = []
gCAIs   = []
for record in CLS_sorted:
    g     = record[0]
    ne    = int(record[1])
    value = CAI_dict.get(g)
    if value == None:
        continue
    gName.append(g)
    nivExpr.append(ne)
    # gCAIs.append(max(value))
    gCAIs.append(np.mean(value))
        
gName   = np.array(gName)
nivExpr = np.array(nivExpr)
gCAIs   = np.array(gCAIs)
#### end strategy 1 ####
'''
# strategy 2 sort nivExpr, note by name of genome, duplicate nivExpr
gName   = []
nivExpr = []
gCAIs   = []

for record in CLS_sorted:
    g     = record[0]
    ne    = int(record[1])
    value = CAI_dict.get(g)
    if value == None:
        continue
    for gCAI in value:
        # gNameNivExprCAIList.append([gName,ne,gCAI])
        gName.append(g)
        nivExpr.append(ne)
        gCAIs.append(gCAI)
        
gName   = np.array(gName)
nivExpr = np.array(nivExpr)
gCAIs   = np.array(gCAIs)
#### end strategy 2 ####

print "Nombre de nivExpr < 30 : ", sum(nivExpr<30)
indice1 = (nivExpr>30)*(nivExpr<100)
print "Nombre de 30 < nivExpr < 100 : ", sum(indice1)
indice2 = nivExpr>100
print "Nombre de nivExpr > 100 : ", sum(indice2)

# saving plots
# print "len(gCAIs): ",gCAIs[:10]
# np.savetxt("gCAIs.nivExpr.txt", np.array(zip(gCAIs,nivExpr)))
plotgCAISNivExpr(gCAIs[indice1],nivExpr[indice1],"30<nivExpr<100")
plotgCAISNivExpr(gCAIs[indice2],nivExpr[indice2],"nivExpr>100")
plotHist(gCAIs,"gCAIsHist")
plotHist(nivExpr[indice1],"30<nivExpr<100Hist")
plotHist(nivExpr[indice2],"nivExpr>100Hist")

'''
# save sorted by gCAIs
M = mergegCAIsNivExprByGname(gCAIs_sorted, CLS_dict)
saveIDgCAInivExpr("gCAIsNivExpr.txt",M)

M  = mergeStep2ByGname(M,step2M)
saveIDgCAInivExprMarked("gCAIsNivExprMark.txt",M)
saveIDgCAInivExprMarked("gCAIsNivExprMarked.txt",M[M[:,3]=='True'])
'''
