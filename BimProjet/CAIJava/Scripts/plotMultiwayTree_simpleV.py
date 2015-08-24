#!env python

import operator
import numpy as np
import pygraphviz as pgv

from BimProjetLib import readClstr, readgCAIs, read2step, getGIDList, domainIdDict, domainGenomeDict, accumulateNivExprDom, sortDictByValue, gCAIsDictToDomDict
from function import GOAnnotation, load_pfam2go, load_goslim, load_godict
from function import plot_switch_view_all_in_one, switch_view, domain_function_list, GOAnnotation, load_all_func_data

#### GLOBAL VARIABLES ####

DF_EDGE_LIST = []
FF_EDGE_LIST = []
FA_EDGE_LIST = []

DOMAIN_SET   = set()
FUNC_SET     = set()
ANCE_SET     = set()

###########################

def reset_global():
    globals()['DF_EDGE_LIST'] = []
    globals()['FF_EDGE_LIST'] = []
    globals()['FA_EDGE_LIST'] = []
    globals()['DOMAIN_SET']   = set()
    globals()['FUNC_SET']     = set()
    globals()['ANCE_SET']     = set()

def sub_load_domain_abundance(domain_list=None):
    # loading
    CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
    # CAI_dict = readgCAIs("../output/cais.lst.2step")
    step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
    gIDList = getGIDList(step2MList)
    # Accumulate Niveau Expression :
    # domIdDict = domainIdDict(gIDList)
    domGenomeDict = domainGenomeDict(gIDList)
    nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
    # nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
    # nivExpr = np.array([ int(val) for val in nivExprAccuDomSortList[:,1]])
    # domains  = nivExprAccuDomSortList[:,0]
    if domain_list <> None:
        ddict = {}
        for d in domain_list:
            ddict[d] = nivExprAccuDomDict[d]
        nivExprAccuDomDict = ddict
    return nivExprAccuDomDict

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

# "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_len30.lst"
# "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_ewvalue_len30.lst"
def sub_load_domain_CAI_dict(fname="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_len30.lst"):
    cais_dict = {}
    CAI_dict = readgCAIs(fname)
    for key,value in CAI_dict.iteritems() :
        key = key.split('_')[-1]
        if cais_dict.has_key(key) :
            cais_dict[key].append(value)
        else:
            cais_dict[key] = [value]
    return cais_dict

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

def find_all_edges(domains, pfam2go_dict, goslimmeta_dict, goDict, search_parent=True):
    for dname in domains:
        if not pfam2go_dict.has_key(dname):
            continue
        go_info_list = pfam2go_dict[dname]
        for go_info in go_info_list:
            if (dname, go_info[2]) not in DF_EDGE_LIST :
                # d = go_info[0]+' '+dname
                d = dname
                DOMAIN_SET.add(d)
                FUNC_SET.add(go_info[2])
                DF_EDGE_LIST.append((d, go_info[2]))
                if search_parent is True:
                    search_parent_and_add(go_info[2], goslimmeta_dict, goDict)
                    
def plot_multiway_tree(plot_name="MultiwayFunctionTree.png", highlight_functions = []):
    G = pgv.AGraph(strict=False,directed=True,ranksep='0.9')
    G.add_edges_from(DF_EDGE_LIST,color='blue')
    G.add_edges_from(FF_EDGE_LIST,color='green')
    G.add_edges_from(FA_EDGE_LIST,color='red')
    G.edge_attr.update(len='13.0')
    for d in DOMAIN_SET :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "cyan"
    for d in ANCE_SET :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "gray"
    for f in highlight_functions :
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "red"
    G.layout()
    G.draw(plot_name)
    return G

def find_most_correlated_functions(seuil):
    f_count_dict = {}
    f_list = []
    for d,f in DF_EDGE_LIST:
        if not f_count_dict.has_key(f):
            f_count_dict[f] = 1
        else:
            f_count_dict[f] += 1

    for key,value in f_count_dict.iteritems():
        if value > seuil:
            f_list.append(key)
            
    return (f_count_dict, f_list)

def filter_domains(f_list):
    d_list = []
    for d,f in DF_EDGE_LIST:
        if f in f_list and not d in d_list :
            d_list.append(d)
    return d_list

def load_all():
    pfam2go_dict = load_pfam2go()
    goslimmeta_dict = load_goslim()
    goDict = load_godict()
    cais_dict = sub_load_domain_CAI_dict()
    domain_abundance = sub_load_domain_abundance(cais_dict.keys())
    return (pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance)
    
def plot_sub_1(seuil=100):
    print "loading... "
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance = load_all()
    domains = np.array(sortDictByValue(domain_abundance))[:,0]
    
    print "searching all edges."
    find_all_edges(domains[:seuil], pfam2go_dict, goslimmeta_dict, goDict)
    
    print "Droping less correlated functions"
    f_count_dict, f_list = find_most_correlated_functions(10)
    
    print "plot..."
    plot_multiway_tree("MultiwayFunctionTreeTOP%dAbondance.png"%seuil, f_list)

def plot_sub_2(seuil=100):    
    print "loading... "
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance = load_all()
    domains = np.array(sortDictByValue(domain_abundance))[:,0]
    
    print "searching all edges."
    find_all_edges(domains, pfam2go_dict, goslimmeta_dict, goDict, False)

    print "Droping less correlated functions"
    f_count_dict, f_list = find_most_correlated_functions(10)

    d_list = filter_domains(f_list)

    print "reset global variable"
    reset_global()
    
    print "searching all edges."
    find_all_edges(d_list[:seuil], pfam2go_dict, goslimmeta_dict, goDict)
    
    print "plot..."
    plot_multiway_tree("MultiwayFunctionTreeTOP%dAbondance.png"%seuil)

def plot_sub_3(seuil=50):    
    print "loading... "
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance = load_all()
    domains = np.array(sortDictByValue(domain_abundance))[:,0]
    
    print "searching all edges."
    find_all_edges(domains[:seuil], pfam2go_dict, goslimmeta_dict, goDict, False)

    print "Droping less correlated functions"
    f_count_dict, f_list = find_most_correlated_functions(10)

    print "filtering domains"
    d_list = filter_domains(f_list)
    
    print "reset global variable"
    reset_global()
    
    print "searching all edges."
    find_all_edges(d_list, pfam2go_dict, goslimmeta_dict, goDict, True)
    
    print "plot..."
    plot_multiway_tree("SimpleMultiwayFunctionTreeTOP%dAbondance.png"%seuil, f_list)

    return d_list
    
def plot_sub_4(seuil=50):    
    print "loading... "
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict_1, domain_abundance = load_all()
    cais_dict_2 = sub_load_domain_CAI_dict("/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_ewvalue_len30.lst")
    domains = np.array(sortDictByValue(domain_abundance))[:,0]

    print "searching all edges."
    find_all_edges(domains[:seuil], pfam2go_dict, goslimmeta_dict, goDict, False)

    print "Droping less correlated functions"
    f_count_dict, f_list = find_most_correlated_functions(10)

    print "filtering domains"
    d_list = filter_domains(f_list)

    dfList = domain_function_list(d_list,pfam2go_dict)
    
    bio_process_1,molec_func_1,cellu_comp_1 = switch_view(dfList,goDict,goslimmeta_dict,cais_dict_1,domain_abundance)
    bio_process_2,molec_func_2,cellu_comp_2 = switch_view(dfList,goDict,goslimmeta_dict,cais_dict_2,domain_abundance)

    plot_switch_view_all_in_one(bio_process_1,bio_process_2,molec_func_1,molec_func_2,cellu_comp_1,cellu_comp_2,"ALLINONETOP%d.png"%seuil,(45,40),caiSeuil=0)

seuil = 400
# plot_sub_1(seuil)
# plot_sub_2(seuil)
# plot_sub_3(seuil)
plot_sub_4(seuil)
