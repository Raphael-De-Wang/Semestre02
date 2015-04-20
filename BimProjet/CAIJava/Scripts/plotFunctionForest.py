#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *
from function import *

import subprocess

def familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict):
    dg = dict()
    de = dict()
    for i,d in enumerate(domains):
        ne = nivExprAccuDomDict.get(d)
        if ne > 500 :
            for g in gCAIsDomDict.get(d):
                if g > 0.8 :
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
CAI_dict = readgCAIs("../output/cais.lst.2step")
step2MList   = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
gIDList = getGIDList(step2MList)

# Accumulate Niveau Expression :
domIdDict = domainIdDict(gIDList)
domGenomeDict = domainGenomeDict(gIDList)
nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
# nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
domains  = nivExprAccuDomSortList[:,0]

# gCAIs par domaine
gCAIsDomDict = gCAIsDictToDomDict(CAI_dict,domIdDict,domains)
# gCAIsDomSortList = np.array(sortDictByValue(CAI_dict))
# save_list('gCAIsDomSortList',gCAIsDomSortList)

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
print len(dfList)
# print len(np.unique(np.array(dfList)[:,1]))

# Construire l'arbre
dot = pydot.Dot(graph_type='digraph',ratio="expand", size='100! 10000!')
for d,func in dfList:
    print "domain : %s, function : %s "%(d,func)
    if func not in goslimmeta_dict.keys():
        add_df_edge(dot,d,func)
    else:
        add_da_edge(dot,d,func)
        
    for key in goDict:
        if func in goDict[key].descendants:
            if key not in goslimmeta_dict.keys():
                # add_ff_edge( dot, func, key)
                search_ancestor(dot,func,key,goDict,goslimmeta_dict)
            else:
                add_fa_edge( dot, func, key)
    
# Dot.write('/tmp/graph.dot', format='raw', prog='dot')
# subprocess.call(["dot", "-Tps", "/tmp/graph.dot", "-o", "/tmp/outfile.ps"])

dot.write_png('/tmp/graph.png', prog='dot')

# print DF_LIST
# print FF_LIST

