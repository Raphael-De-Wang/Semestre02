#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *
from function import *

import subprocess

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
    
# loading
CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
CAI_dict = readgCAIs("../output/cais.lst.2step")
step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
gIDList = getGIDList(step2MList)

# mean,std of CAI
caiList = []
for key,CAI in CAI_dict.iteritems():
    caiList.append(CAI)
moyCAI = np.mean(caiList)
stdCAI = np.std(caiList)

# Accumulate Niveau Expression :
domIdDict = domainIdDict(gIDList)
domGenomeDict = domainGenomeDict(gIDList)
nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
# nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
domains  = nivExprAccuDomSortList[:,0]

# mean,std of Niveau Expression :
neList = []
for key,ne in nivExprAccuDomDict.iteritems():
    neList.append(ne)
moyNE = np.mean(neList)
stdNE = np.std(neList)
    
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

# deSortList,dg,de = familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict,neSeuil=moyNE,gSeuil=moyCAI)
deSortList,dg,de = familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict,neSeuil=0,gSeuil=0)
print 'Domain Expression Level List Done.'
#### 
klist = []
vseuil = 0.6
ndg = {}
nde = {}
for k,vl in dg.iteritems():
    if np.max(vl) > vseuil:
        klist.append(k)
for k in klist:
    ndg[k] = dg[k]
    nde[k] = de[k]
deSortList = np.array(sortDictByValue(nde))
print 'List Filtering Done. '
####
# domainHTML(deSortList,pfam2go_dict,dg,de)
dfList = domain_function_list(np.array(deSortList)[:,0],pfam2go_dict)
# dfDict = domain_function_dict(np.array(deSortList)[:,0],pfam2go_dict)
# print len(dfList)
# print len(np.unique(np.array(dfList)[:,1]))
print 'Domain to GoTerms Done. '
print 'Domain Function List Length: ', len(dfList)

# dfList = [['PF01121', 'GO:0015937']]

def constr_arbre():
    # Construire l'arbre
    dot = pydot.Dot(graph_type='digraph' ,ratio="expand", size='100! 10000!')
    for d,func in reversed(dfList):
        # print "domain : %s, function : %s "%(d,func)
        if func not in goslimmeta_dict.keys():
            add_df_edge(dot,d,func)
        else:
            add_da_edge(dot,d,func)
            
        for key in goDict:
            if func in goDict[key].descendants:
                if key not in goslimmeta_dict.keys():
                    add_ff_edge( dot, func, key) ####
                    search_ancestor(dot,func,key,goDict,goslimmeta_dict)
                else:
                    add_fa_edge( dot, func, key)
    
    # Dot.write('/tmp/graph.dot', format='raw', prog='dot')
    # subprocess.call(["dot", "-Tps", "/tmp/graph.dot", "-o", "/tmp/outfile.ps"])
    print 'writing plot file. '
    dot.write_png('/tmp/graph.png', prog='dot')

print 'Building Function tree. '
# constr_arbre()
print "length of DF_LIST : ", len(DF_LIST)
print "length of FF_LIST : ", len(FF_LIST)
print "DF_LIST : ", DF_LIST
print "FF_LIST : ", FF_LIST

def plot_nbr_domain_par_func(dfDict,width=0.35,fname=None):
    dfList = np.array([ [k,v] for k,v in dfDict.iteritems() ])
    ind = np.arange(len(dfList))
    fig = plt.figure(figsize=(100,7))
    fig.subplots_adjust(bottom=0.3)
    plt.title('Nombre de domain par function, CAI > %f, Niveau D\'Expression > %f'%(moyCAI,moyNE))
    funcList = dfList[:,0]
    nbrList  = [ int(n) for n in dfList[:,1] ]
    plt.bar(ind, nbrList, width, color="#2aa198")
    plt.xticks(ind+width/4.,funcList,rotation=80)
    if fname==None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)

# plot_nbr_domain_par_func(dfDict,fname='nbr_domain_par_func.png')

# bio_process,molec_func,cellu_comp = switch_view(dfList,goDict,goslimmeta_dict,dg,de)
bio_process,molec_func,cellu_comp = switch_view(dfList,goDict,goslimmeta_dict,ndg,nde)

# plot_switch_view('Biological Process',bio_process,'BiologicalProcess.png',(15,50))
# plot_switch_view('Molecular Function',molec_func, 'MolecularFunction.png',(45,30))
# plot_switch_view('Cellular Component',cellu_comp, 'CellularComponent.png',(15,15))

plot_switch_view('Biological Process',bio_process,'BiologicalProcess_meta.png',(30,15))
plot_switch_view('Molecular Function',molec_func, 'MolecularFunction_meta.png',(40,25))
plot_switch_view('Cellular Component',cellu_comp, 'CellularComponent_meta.png',(30,15))

