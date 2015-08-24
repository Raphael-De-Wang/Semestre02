#!env python

import numpy as np
import pygraphviz as pgv

from BimProjetLib import readClstr, readgCAIs, read2step, getGIDList, domainIdDict, domainGenomeDict, accumulateNivExprDom, sortDictByValue, gCAIsDictToDomDict
from function import GOAnnotation, load_pfam2go, load_goslim, load_godict

DF_EDGE_LIST = []
FF_EDGE_LIST = []
FA_EDGE_LIST = []

DOMAIN_SET   = set()
FUNC_SET     = set()
ANCE_SET     = set()

def sub_load_sorted_by_abundance():
    # loading
    CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
    # CAI_dict = readgCAIs("../output/cais.lst.2step")
    step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
    gIDList = getGIDList(step2MList)
    # Accumulate Niveau Expression :
    # domIdDict = domainIdDict(gIDList)
    domGenomeDict = domainGenomeDict(gIDList)
    nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
    nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
    # nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
    # domains  = nivExprAccuDomSortList[:,0]
    return nivExprAccuDomSortList

def sub_load_sorted_by_CAI():
    cais_max_dict = {}
    # CAI_dict = readgCAIs("../output/cais.lst.2step")
    CAI_dict = readgCAIs("../output/cais_len30.lst")
    for key,value in CAI_dict.iteritems() :
        key = key.split('_')[-1]
        if cais_max_dict.has_key(key) and cais_max_dict[key] > value :
            continue
        else:
            cais_max_dict[key] = value
    return np.array(sortDictByValue(cais_max_dict))

def is_ancestor(goterm, goslimmeta_dict) :
    if goterm in goslimmeta_dict.keys() :
        return True
    return False

def is_descendant(goterm, go_dict, key):
    if goterm in go_dict[key].descendants:
        return True
    return False

def search_parent_and_add(goterm, goslimmeta_dict, goDict):
    for k,v in goDict.iteritems():
        if is_descendant(goterm, goDict, k):
            if [goterm, k] not in FF_EDGE_LIST + FA_EDGE_LIST:
                if not is_ancestor(k, goslimmeta_dict) :
                    FF_EDGE_LIST.append([goterm,k])
                    FUNC_SET.add(k)
                    search_parent_and_add(k, goslimmeta_dict, goDict)
                else:
                    FA_EDGE_LIST.append([goterm,k])
                    ANCE_SET.add(k)

def find_all_edges(domains, pfam2go_dict, goslimmeta_dict, goDict):
    for dname in domains:
        if not pfam2go_dict.has_key(dname):
            continue
        go_info_list = pfam2go_dict[dname]
        for go_info in go_info_list:
            if (dname, go_info[2]) not in DF_EDGE_LIST :
                d = go_info[0]+' '+dname
                DOMAIN_SET.add(d)
                FUNC_SET.add(go_info[2])
                DF_EDGE_LIST.append((d, go_info[2]))
                search_parent_and_add(go_info[2], goslimmeta_dict, goDict)
                    
def plot_multiway_tree(plot_name="MultiwayFunctionTree.png"):
    G = pgv.AGraph(strict=False,directed=True,ranksep='0.9')
    G.add_edges_from(DF_EDGE_LIST,color='blue')
    G.add_edges_from(FF_EDGE_LIST,color='green')
    G.add_edges_from(FA_EDGE_LIST,color='red')
    G.edge_attr.update(len='20.0')
    for d in DOMAIN_SET :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "cyan"
    G.layout()
    G.draw(plot_name)
    return G

pfam2go_dict = load_pfam2go()
goslimmeta_dict = load_goslim()
goDict = load_godict()

domains = sub_load_sorted_by_abundance()[:,0]
# domains = sub_load_sorted_by_CAI()[:,0]

find_all_edges(domains[:100], pfam2go_dict, goslimmeta_dict, goDict)
plot_multiway_tree("MultiwayFunctionTreeTOP100Abondance.png")

