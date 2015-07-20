#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from BimProjetLib import *
from function import *

import subprocess

def familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict,neSeuil=500,gSeuil=0.8):
    dg = dict()
    de = dict()
    for i,d in enumerate(domains):
        ne = nivExprAccuDomDict.get(d)
        if ne > neSeuil :
            for g in gCAIsDomDict.get(d):
                if g > gSeuil :
                    if dg.has_key(d):
                        dg[d].append(g)
                    else:
                        dg[d] = [g]
                    if not de.has_key(d):
                        de[d] = ne
                    elif de[d] <> ne :
                            raise ValueError("Expression Level Value Conflict.")
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
    
def domain_function_dict(domains,pfam2go_dict):
    dfDict = {}
    for i,d in enumerate(domains):
        if pfam2go_dict.has_key(d):
            for record in pfam2go_dict[d]:
                func = record[2]
                if dfDict.has_key(func):
                    dfDict[func] += 1
                else:
                    dfDict[func] = 1
    return dfDict
    
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
vseuil = 0.9
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

def constr_arbre():
    # Construire l'arbre
    dot = pydot.Dot(graph_type='digraph' ,ratio="expand", size='100! 10000!')
    for d,func in dfList:
        print "domain : %s, function : %s "%(d,func)
        if len(DOMAIN_LIST)%100 == 0:
            print DOMAIN_LIST
        if len(FUNC_LIST)%100 == 0:
            print FUNC_LIST
        if len(DF_LIST)%100 == 0:
            print DF_LIST
        if len(FF_LIST)%100 == 0:
            print FF_LIST
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
constr_arbre()
# print DF_LIST
# print FF_LIST
exit()

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

def add_in_dict(D,key,name,dom,cais,ne):
    if not D.has_key(key):
        D[key] = (name,[dom],[cais],[ne])
    else:
        n,dom_list,cais_list,ne_list = D[key]
        if n <> name:
            raise ValueError('different function name : %s, %s'%(n,name))
        dom_list.append(dom)
        cais_list.append(cais)
        ne_list.append(ne)
        D[key] = (name,dom_list,cais_list,ne_list)

def switch_view(dfList,goDict,goslimmeta_dict,dgDict,deDict):
    bio_process = {}
    molec_func  = {}
    cellu_comp  = {}
    for d,func in dfList:
        if func in goslimmeta_dict.keys():
            record = goslimmeta_dict.get(func)
        elif func in goDict.keys():
            record = goDict.get(func)
        else:
            continue
        ne  = deDict[d]
        cais= dgDict[d]
        name= record.name
        if record.type == 'biological_process':
            add_in_dict(bio_process,func,name,d,cais,ne)
        if record.type == 'molecular_function': 
            add_in_dict(molec_func,func,name,d,cais,ne)
        if record.type == 'cellular_component': 
            add_in_dict(cellu_comp,func,name,d,cais,ne)
    return bio_process,molec_func,cellu_comp
    
# bio_process,molec_func,cellu_comp = switch_view(dfList,goDict,goslimmeta_dict,dg,de)
bio_process,molec_func,cellu_comp = switch_view(dfList,goDict,goslimmeta_dict,ndg,nde)

def plot_switch_view(func_type,func_dict,fname=None,figsize=None,caiSeuil=0.0):
    nameList = []
    domList  = []
    caisMinList = []
    caisMaxList = []
    caisMeanList= []
    neList      = []
    funcList    = []
    fndict      = {}
    for key,value in func_dict.iteritems():
        (name,dom_list,cais_list,ne_list) = value
        if not fndict.has_key(key):
            fndict[key] = sum(ne_list)
        else :
            fndict[key] += sum(ne_list)
    for key,mNe in sortDictByValue(fndict,False):
        value = func_dict[key]
        (name,dom_list,cais_list,ne_list) = value
        cais_list = [j for i in cais_list for j in i]
        if max(cais_list) < caiSeuil :
            continue
        caisMinList.append(min(cais_list))
        caisMaxList.append(max(cais_list))
        caisMeanList.append(np.mean(cais_list))
        nameList.append(name)
        domList.append(len(dom_list))
        neList.append(sum(ne_list))
        funcList.append(key)
    bottoms = np.arange(len(nameList))*1.2
    fig = plt.figure(figsize=figsize)
    plt.suptitle(func_type, fontsize=15)
    ax = plt.subplot(131,axisbg="#fdf6e3")
    ax.set_title("Nombre de Domain par Function")
    plt.bar(left=np.zeros(len(nameList)),
            width=domList, bottom=bottoms,align="center",
            color="#2aa198",orientation="horizontal",height=1.0)
    plt.yticks(bottoms+0.1,nameList)
    ax = plt.subplot(132,axisbg="#fdf6e3")
    ax.set_title("Niveau d'Expression")
    plt.bar(left=np.zeros(len(nameList)),
            width=neList, bottom=bottoms,align="center",
            color="#2aa198",orientation="horizontal",height=1.0)
    plt.yticks(bottoms+0.1,funcList)
    ax = plt.subplot(133,axisbg="#fdf6e3")
    ax.set_title("gCAI")
    plt.bar(left=np.zeros(len(nameList)),align="center",
            width=caisMaxList, bottom=bottoms,
            color="b",orientation="horizontal",
            height=0.8,label='Max')
    plt.bar(left=np.zeros(len(nameList)),align="center",
            width=caisMeanList, bottom=bottoms,
            color="g",orientation="horizontal",
            height=0.9,label='Mean')
    plt.bar(left=np.zeros(len(nameList)),align="center",
            width=caisMinList, bottom=bottoms,
            color="r",orientation="horizontal",
            height=1.0,label='Min')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.yticks(bottoms+0.1,funcList)
    # ax.yaxis.tick_right()
    # fig.subplots_adjust(left=0.3,hspace=0.05,bottom=0.1,top=0.90)
    fig.subplots_adjust(left=0.3,bottom=0.1,top=0.90)
    if fname == None:
        plt.show()
    else:
        plt.savefig(fname)
        handle = open(fname.replace('png','txt'),'w')
        for f in funcList:
            handle.write("%s\n"%f)
        handle.close()
    plt.close(fig)

# plot_switch_view('Biological Process',bio_process,'BiologicalProcess.png',(15,50))
# plot_switch_view('Molecular Function',molec_func, 'MolecularFunction.png',(45,30))
# plot_switch_view('Cellular Component',cellu_comp, 'CellularComponent.png',(15,15))

plot_switch_view('Biological Process',bio_process,'BiologicalProcess.png',(30,15))
plot_switch_view('Molecular Function',molec_func, 'MolecularFunction.png',(40,25))
plot_switch_view('Cellular Component',cellu_comp, 'CellularComponent.png',(30,15))
