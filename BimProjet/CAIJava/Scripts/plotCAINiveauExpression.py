#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *

# strategy 1 mean ou max of duplications, func: mean(), max(), ...
def sortByClassJoinCAIs(sortCLS,CAI_dict,func):
    gName   = []
    nivExpr = []
    gCAIs   = []
    for record in sortCLS:
        g     = record[0]
        ne    = int(record[1])
        value = CAI_dict.get(g)
        if value == None:
            continue
            gName.append(g)
            nivExpr.append(ne)
            gCAIs.append(func(value))
    gName   = np.array(gName)
    nivExpr = np.array(nivExpr)
    gCAIs   = np.array(gCAIs)
    return (gName,nivExpr,gCAIs)
#### end strategy 1 ####

# strategy 2 sort nivExpr, note by name of genome, duplicate nivExpr
def sortByClassDupCAIs(sortCLS,CAI_dict):
    gName   = []
    nivExpr = []
    gCAIs   = []
    for record in sortCLS:
        g     = record[0]
        ne    = int(record[1])
        value = CAI_dict.get(g)
        if value == None:
            continue
        for gCAI in value:
            gName.append(g)
            nivExpr.append(ne)
            gCAIs.append(gCAI)
    gName   = np.array(gName)
    nivExpr = np.array(nivExpr)
    gCAIs   = np.array(gCAIs)
    return (gName,nivExpr,gCAIs)
#### end strategy 2 ####

#### strategy 3 ####
def sortByCAIsJoinNivExprDom(sortCAIs,nivExprDomDict):
    gNames  = []
    nivExpr = []
    gCAIs   = []
    for record in sortCAIs:
        g     = record[0]
        dn    = g.split('_')[-1]
        cai   = float(record[1])
        ne    = nivExprDomDict.get(dn)
        if ne == None:
            continue
        nivExpr.append(ne)
        gCAIs.append(cai)
        gNames.append(g)
    nivExpr = np.array(nivExpr)
    gCAIs   = np.array(gCAIs)
    gNames  = np.array(gNames)
    return (gNames,nivExpr,gCAIs)
    
#### end strategy 3 ####

#### strategy 4 ####
def sortByNivExprDomJoinCAIs(sortCAIs,nivExprDomDict):
    gNames,nivExpr,gCAIs = sortByCAIsJoinNivExprDom(sortCAIs,nivExprDomDict)
    lstSort = sortListByCol(zip(gNames,nivExpr,gCAIs),1)
    gNames,nivExpr,gCAIs = [],[],[]
    for record in lstSort:
        gNames.append(record[0])
        nivExpr.append(record[1])
        gCAIs.append(record[2])
    return np.array(gNames),np.array(nivExpr),np.array(gCAIs)
#### end strategy 4 ####

# loading
CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
CAI_dict = readgCAIs("../output/cais.lst.bestdomain.e-8")
step2MList   = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
nivExprDomDict = nivExprDomain(step2MList, 4)
# step2MList[binSearch(step2MList,step2GnameComp,"GG7SD3401DYJ3N")]

nivExprDomDict = nivExprDomain(step2MList,4)
print "Domaine Nombre: ",len(nivExprDomDict)
# sortCLS = np.array(sortDictByValue(CLS_dict))
sortCAIs = np.array(sortDictByValue(CAI_dict))

# gID,nivExpr,gCAIs = sortByCAIsJoinNivExprDom(sortCAIs,nivExprDomDict)
gID,nivExpr,gCAIs = sortByNivExprDomJoinCAIs(sortCAIs,nivExprDomDict)

# Accumulate Niveau Expression :
domIdDict = domainIdDict(gID)
nivExprAccuDomDict = accumulateNivExprDom(nivExprDomDict,domIdDict)
nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
domains  = nivExprAccuDomSortList[:,0]
plotBarChart(nivExpr[:100],domains[:100])

# gCAIs par domaine
gCAIsDomDict = gCAIsDictToDomDict(CAI_dict,domIdDict,domains)
gCAIsDomList = np.array([[i,g] for i,d in enumerate(domains) for g in gCAIsDomDict.get(d) ])
plotgCAIsDomain(gCAIsDomList[:20000,1],gCAIsDomList[:20000,0])

print "Nombre de nivExpr < 30 : ", sum(nivExpr<30)
indice1 = (nivExpr>30)*(nivExpr<100)
print "Nombre de 30 < nivExpr < 100 : ", sum(indice1)
indice2 = nivExpr>100
print "Nombre de nivExpr > 100 : ", sum(indice2)


# saving plots
# np.savetxt("gCAIs.nivExpr.txt", np.array(zip(gCAIs,nivExpr)))
plotgCAISNivExpr(gCAIs[indice1],nivExpr[indice1],"30<nivExpr<100")
plotgCAISNivExpr(gCAIs[indice2],nivExpr[indice2],"nivExpr>100")
plotHist(gCAIs,"gCAIsHist")
'''
plotHist(nivExpr[nivExpr<30],"nivExpr<30Hist")
plotHist(nivExpr[indice1],"30<nivExpr<100Hist")
plotHist(nivExpr[indice2],"nivExpr>100Hist")
'''
plotHist(nivExpr,"nivExprHist")

# save sorted by gCAIs
# M = mergegCAIsNivExprByGname(gCAIs_sorted, CLS_dict)
# saveIDgCAInivExpr("resultats/gCAIsNivExpr.txt",np.array(zip(gID,gCAIs,nivExpr)))

'''
M  = mergeStep2ByGname(M,step2MList)
saveIDgCAInivExprMarked("gCAIsNivExprMark.txt",M)
saveIDgCAInivExprMarked("gCAIsNivExprMarked.txt",M[M[:,3]=='True'])
'''

