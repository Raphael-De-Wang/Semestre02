#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *
from function import *

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

def familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict):
    dg = dict()
    de = dict()
    for i,d in enumerate(domains):
        ne = nivExprAccuDomDict.get(d)
        if ne > 0 :
            for g in gCAIsDomDict.get(d):
                if g > 0.6 :
                    if dg.has_key(d):
                        dg[d].append(g)
                    else:
                        dg[d] = [g]
                    if not de.has_key(d):
                        de[d] = ne
    deSortList = np.array(sortDictByValue(de))
    return deSortList,dg,de

def domainHTML(deSortList,pfam2go_dict,dg,de):
    print '<!DOCTYPE html>'
    print '<html>'
    print '<body>'
    print '<table style="width:40%">'
    print '<tr><td>Domain</td><td>Expression Level</td><td>gCAI</td></tr>'
    for d,e in deSortList:
        dg[d].sort(reverse=True)
        print '<tr>'
        print '<td><a href="http://pfam.xfam.org/family/%s">'%d,d,'</a></td>'
        print '<td>',e,'</td>'
        print '<td>',dg[d][0],'</td>'
        if pfam2go_dict.has_key(d):
            for record in pfam2go_dict[d]:
                func = record[2]
                print '<td><a href="http://www.ebi.ac.uk/QuickGO/GTerm?id=%s">'%func,func,'</a></td>'
        print '</tr>'
    print '</table>'
    print '</body>'
    print '</html>'
    
def domain_function_list(domains,pfam2go_dict):
    dfList = []
    for i,d in enumerate(domains):
        if pfam2go_dict.has_key(d):
            for record in pfam2go_dict[d]:
                func = record[2]
                dfList.append([d,func])
    return dfList
    
# loading
CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
# CAI_dict = readgCAIs("../output/cais.lst.2step")
CAI_dict = readgCAIs("../output/AT_arc_metatrans.filtered.fasta.cais_list")
step2MList   = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
# step2MList[binSearch(step2MList,step2GnameComp,"GG7SD3401DYJ3N")]

# nivExprDomDict = nivExprDomain(step2MList,4)
# print "Domaine Nombre: ",len(nivExprDomDict)
# sortCLS = np.array(sortDictByValue(CLS_dict))
# sortCAIs = np.array(sortDictByValue(CAI_dict))

# gID,nivExpr,gCAIs = sortByCAIsJoinNivExprDom(sortCAIs,nivExprDomDict)
# gID,nivExpr,gCAIs = sortByNivExprDomJoinCAIs(sortCAIs,nivExprDomDict)
gIDList = getGIDList(step2MList)

# Accumulate Niveau Expression :
domIdDict = domainIdDict(gIDList)
domGenomeDict = domainGenomeDict(gIDList)
nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
domains  = nivExprAccuDomSortList[:,0]

# save_list('nivExprAccuDomSortList', nivExprAccuDomSortList)

'''
moy = np.mean(nivExpr)
std = np.std(nivExpr)
plotBarChart(nivExpr,domains,"NivExprAll")
plotBarChart(nivExpr[nivExpr<moy],domains[nivExpr<moy],"nivExpr<moy")
plotBarChart(nivExpr[nivExpr<moy+std],domains[nivExpr<moy+std],"nivExprmoy+std")
plotBarChart(nivExpr[nivExpr<moy+2*std],domains[nivExpr<moy+2*std],"nivExprmoy+2std")
plotBarChart(nivExpr[nivExpr<moy+3*std],domains[nivExpr<moy+3*std],"nivExprmoy+3std")
plotBarChart(nivExpr[nivExpr<moy+4*std],domains[nivExpr<moy+4*std],"nivExprmoy+4std")
'''

# gCAIs par domaine
gCAIsDomDict = gCAIsDictToDomDict(CAI_dict,domIdDict,domains)
gCAIsDomSortList = np.array(sortDictByValue(CAI_dict))
# save_list('gCAIsDomSortList',gCAIsDomSortList)

xDiv=1
yDiv=20

gCAIsDomList = np.array([[i,g] for i,d in enumerate(domains) for g in gCAIsDomDict.get(d) if type(g) == float ])
# indice = gCAIsDomList[:,0] < top
plotgCAIsDomain(gCAIsDomList[:,1],gCAIsDomList[:,0],yDiv,xDiv,figName="gCAIsChaleur[%d][xDiv=%d][yDiv=%d]"%(len(np.unique(gCAIsDomList[:,0])),xDiv,yDiv))
# exit()

plotgCAIsDomain2(gCAIsDomList[:,0],gCAIsDomList[:,1],"gCAIsDistr")
exit()
'''
indice = gCAIsDomList[:,0] < 1000
plotgCAIsDomain(gCAIsDomList[indice,1],gCAIsDomList[indice,0],"gCAIsChaleurTop1000")
indice = gCAIsDomList[:,0] < 500
plotgCAIsDomain(gCAIsDomList[indice,1],gCAIsDomList[indice,0],"gCAIsChaleurTop500")
plotgCAIsDomain2(gCAIsDomList[:,1],gCAIsDomList[:,0],"gCAIsDistr")

gCAIsDomNivExprList = np.array([[i,g,nivExprAccuDomDict.get(d)] for i,d in enumerate(domains) for g in gCAIsDomDict.get(d) ])
plotgCAIsDomain3(gCAIsDomNivExprList[:,1],gCAIsDomNivExprList[:,0],gCAIsDomNivExprList[:,2],"gCAIsDistrColoreExpr")
'''
# gCAIsDomNivExprList = np.array([[g,nivExprAccuDomDict.get(d),d] for i,d in enumerate(domains) for g in gCAIsDomDict.get(d) if g > 0.8 and nivExprAccuDomDict.get(d) > 1000 ])
# plotgCAIsDomain4(gCAIsDomNivExprList[:,1],gCAIsDomNivExprList[:,0],xDiv,yDiv,figName="gCAIsExprChaleur[xDiv=%d][yDiv=%d]"%(xDiv,yDiv))
# print np.array([ [d,nivExprAccuDomDict.get(d),g] for i,d in enumerate(domains) for g in gCAIsDomDict.get(d) if g > 0.8 and nivExprAccuDomDict.get(d) > 1000 ])

handle = open("pfam2go")
pfam2go_dict = load_pfam2go_toDict(handle)
handle.close()

# fd_dict = pfam2go_to_funcDomainDict(pfam2go_dict)

handle = open("GO_SLIM_META_DICT.txt",'r')
goslimmeta_dict = loadToDict(handle)
handle.close()

handle = open("GO_DICT.txt",'r')
goDict = loadToDict(handle)
handle.close()

deSortList,dg,de = familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict)
# domainHTML(deSortList,pfam2go_dict,dg,de)
dfList = domain_function_list(np.array(deSortList)[:,0],pfam2go_dict)
# print len(np.array(dfList))
# print len(np.unique(np.array(dfList)[:,1]))

step = 100
for i in range(100,401,step):
    indice = ( gCAIsDomList[:,0] < i+1)*(gCAIsDomList[:,0] > i - step )
    plotgCAIsDomain(gCAIsDomList[indice,1],gCAIsDomList[indice,0],"gCAIsChaleur%d"%i)
    plotgCAIsDomain2(gCAIsDomList[indice,1],gCAIsDomList[indice,0]-i+step+1,"gCAIsDistr%d"%i)
exit()
print "Nombre de nivExpr < 30 : ", sum(nivExpr<30)
indice1 = (nivExpr>30)*(nivExpr<100)
print "Nombre de 30 < nivExpr < 100 : ", sum(indice1)
indice2 = nivExpr>100
print "Nombre de nivExpr > 100 : ", sum(indice2)

exit()
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

