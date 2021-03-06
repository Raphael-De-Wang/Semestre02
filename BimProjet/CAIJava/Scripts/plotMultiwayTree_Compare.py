#!env python

import csv
import cPickle
import operator
import numpy as np
import pygraphviz as pgv

from BimProjetLib import readClstr, readgCAIs, read2step, getGIDList, domainIdDict, domainGenomeDict, accumulateNivExprDom, sortDictByValue, gCAIsDictToDomDict
from function import GOAnnotation, load_pfam2go, load_goslim, load_godict
from function import plot_switch_view_all_in_one, switch_view, domain_function_list, GOAnnotation, load_all_func_data, load_domain_gcai_abundance

#### GLOBAL VARIABLES ####

DF_EDGE_LIST = []
FF_EDGE_LIST = []
FA_EDGE_LIST = []

DOMAIN_SET   = set()
FUNC_SET     = set()
ANCE_SET     = set()

###########################

def reset_global_list():
    globals()['DF_EDGE_LIST'] = []
    globals()['FF_EDGE_LIST'] = []
    globals()['FA_EDGE_LIST'] = []

def reset_global_set():    
    globals()['DOMAIN_SET']   = set()
    globals()['FUNC_SET']     = set()
    globals()['ANCE_SET']     = set()
    
def reset_global():
    reset_global_list()
    reset_global_set()
    
def copy_global_list():
    DF_EDGE_LIST = globals()['DF_EDGE_LIST']
    FF_EDGE_LIST = globals()['FF_EDGE_LIST'] 
    FA_EDGE_LIST = globals()['FA_EDGE_LIST']
    return (DF_EDGE_LIST, FF_EDGE_LIST, FA_EDGE_LIST)
    
def copy_global_set():
    DOMAIN_SET   = globals()['DOMAIN_SET']
    FUNC_SET     = globals()['FUNC_SET']   
    ANCE_SET     = globals()['ANCE_SET']
    return (DOMAIN_SET, FUNC_SET, ANCE_SET)

def sub_load_domain_abundance(domain_list=None):
    # loading
    CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
    step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
    gIDList = getGIDList(step2MList)
    domGenomeDict = domainGenomeDict(gIDList)
    nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
    if domain_list <> None:
        ddict = {}
        for d in domain_list:
            ddict[d] = nivExprAccuDomDict[d]
        nivExprAccuDomDict = ddict
    return nivExprAccuDomDict

def sub_load_sorted_by_CAI():
    cais_max_dict = {}
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
                    search_parent_and_add(k, goslimmeta_dict, goDict)
                else:
                    FA_EDGE_LIST.append([goterm,k])
    for f1, f2 in FF_EDGE_LIST :
        FUNC_SET.add(f1)
        FUNC_SET.add(f2)
    for f,a in FA_EDGE_LIST :
        FUNC_SET.add(f)
        ANCE_SET.add(a)
        
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
                    
def find_all_fonction_edges(fonctions, goslimmeta_dict, goDict):
    for f in fonctions:
        FUNC_SET.add(f)
        search_parent_and_add(f, goslimmeta_dict, goDict)
                    
def plot_multiway_tree(plot_name, df_edge_list, ff_edge_list, fa_edge_list,
                       domain_set_1, func_set_1, ance_set_1,
                       domain_set_2=set(), func_set_2=set(), ance_set_2=set(),
                       func_commun=[], ance_commun=[], base_commun=[],
                       children_1=[], children_2=[],
                       label_func_list=[], label_dict = {}):
    
    G = pgv.AGraph(strict=False,directed=True,ranksep='0.9')
    
    G.add_edges_from(df_edge_list,color='#9e9eff') # blue
    G.add_edges_from(ff_edge_list,color='#9dff9d') # green
    G.add_edges_from(fa_edge_list,color='#ff9d9d') # red
    G.edge_attr.update(len='20.0')

    for d in domain_set_1 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "cyan"
    for d in func_set_1 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "lightblue"
    for d in ance_set_1 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "gray"
        
    for d in domain_set_2 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "cyan"
    for d in func_set_2 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "#cbffb1" # "palegreen"
    for d in ance_set_2 :
        n = G.get_node(d)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "gray"

    for f in children_1:
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "blue"
    for f in children_2 :
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "#14ff14" # "green"

    for f in func_commun :
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "pink"
    for f in ance_commun :
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "orange"

    for f in base_commun :
        n = G.get_node(f)
        n.attr["style"] = "filled"
        n.attr["fillcolor"] = "red"

    for f in label_func_list :
        n = G.get_node(f)        
        n.attr["label"] = label_dict[f]
        
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

def load_all_metagenomic():
    pfam2go_dict = load_pfam2go()
    goslimmeta_dict = load_goslim()
    goDict = load_godict()
    cais_dict = sub_load_domain_CAI_dict()
    domain_abundance = sub_load_domain_abundance(cais_dict.keys())
    return (pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance)

def load_all_opposite():
    # opposite organismes, gcais are calculate by signature metagenomic
    fname = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline/domain_gcai_abundance_oppo_by_oppo_moy.csv"
    # pfam2go_dict, goslimmeta_dict, goDict, dom_gcai, dom_abundance  = load_all_func_data()
    dom_gcai_oppo, dom_abundance_oppo = load_domain_gcai_abundance(fname)
    return dom_gcai_oppo, dom_abundance_oppo

def merge_dict(d1,d2):
    dm = d2
    for k,v in d1.iteritems():
        if not dm.has_key(k):
            dm[k] = v
        else :
            dm[k] = dm[k] + d1[k]
    return dm
    
def load_all_meta_by_meta():
    pfam2go_dict = load_pfam2go()
    goslimmeta_dict = load_goslim()
    goDict = load_godict()
    cais_dict = sub_load_domain_CAI_dict()
    domain_abundance = sub_load_domain_abundance(cais_dict.keys())
    fname = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline/domain_gcai_abundance_meta_all_by_meta.csv"
    dom_gcai_meta, dom_abundance_meta = load_domain_gcai_abundance(fname)
    dom_gcai = merge_dict(cais_dict,dom_gcai_meta)
    dom_abundance = merge_dict(domain_abundance,dom_abundance_meta)
    return (pfam2go_dict, goslimmeta_dict, goDict, dom_gcai, dom_abundance)
    
def load_all_meta_by_meta_moy():
    pfam2go_dict = load_pfam2go()
    goslimmeta_dict = load_goslim()
    goDict = load_godict()
    cais_dict = sub_load_domain_CAI_dict("/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/AT_arc_metatrans_gcai_meta_moy.csv")
    domain_abundance = sub_load_domain_abundance(cais_dict.keys())
    fname = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline/domain_gcai_abundance_meta_all_by_meta_moy.csv"
    dom_gcai_meta, dom_abundance_meta = load_domain_gcai_abundance(fname)
    dom_gcai = merge_dict(cais_dict,dom_gcai_meta)
    dom_abundance = merge_dict(domain_abundance,dom_abundance_meta)
    return (pfam2go_dict, goslimmeta_dict, goDict, dom_gcai, dom_abundance)
    
def fonction_abondance(pfam2go_dict, domain_abundance):
    go_ab_dict = {}
    for d,ab in domain_abundance.iteritems() :
        if pfam2go_dict.has_key(d):
            for func_name, go_name, go_id in pfam2go_dict[d]:
                if go_ab_dict.has_key(go_id) :
                   go_ab_dict[go_id] += int(ab)
                else:
                   go_ab_dict[go_id] = int(ab)
    return go_ab_dict
    
def fonction_max_gcai(pfam2go_dict, domain_cais_dict):
    go_gcai_dict = {}
    for d,cais in domain_cais_dict.iteritems() :
        if pfam2go_dict.has_key(d):
            for func_name, go_name, go_id in pfam2go_dict[d]:
                if not go_gcai_dict.has_key(go_id) or go_gcai_dict[go_id] < max(cais) :
                    go_gcai_dict[go_id] = max(cais)
    return go_gcai_dict

def filter_go_ab_dict_by_gcai(go_ab_dict, go_gcai_dict, gcai_threshold):
    go_ab_dict_clean = {}
    for go,ab in go_ab_dict.iteritems():
        if go_gcai_dict[go] >= gcai_threshold :
            go_ab_dict_clean[go] = ab
    return go_ab_dict_clean
    
def plot_sub_1(range1=(0,50), range2=(50,100)):
    print "loading..."
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict, domain_abundance = load_all_metagenomic()
    
    print "convert to fonction domain list..."
    go_ab_dict = fonction_abondance(pfam2go_dict, domain_abundance)
    go_gcai_dict = fonction_max_gcai(pfam2go_dict, cais_dict)

    print "filtering by gcai threshold..."
    gcai_threshold = 0.6
    go_ab_dict_clean = filter_go_ab_dict_by_gcai(go_ab_dict, go_gcai_dict, gcai_threshold)

    if len(go_ab_dict_clean.keys()) < max(range1+range2):
        print "Number of funtion is less than max of range requist."
        exit()
    
    print "sort functions by abondance..."
    go_ab_list = sortDictByValue(go_ab_dict_clean)
    func_list  = np.array(go_ab_list)[:,0]

    print "find all parents of range1..."
    reset_global()
    find_all_fonction_edges(func_list[range1[0]:range1[1]], goslimmeta_dict, goDict)
    DOMAIN_SET_1, FUNC_SET_1, ANCE_SET_1 = copy_global_set()
    
    print "find all parents of range2..."
    reset_global()
    find_all_fonction_edges(func_list[range2[0]:range2[1]], goslimmeta_dict, goDict)
    DOMAIN_SET_2, FUNC_SET_2, ANCE_SET_2 = copy_global_set()
    
    print "searching commun functions between range 1 and range 2..."
    func_commun = [ f for f in FUNC_SET_1 if f in FUNC_SET_2 ]
    ance_commun = [ f for f in ANCE_SET_1 if f in ANCE_SET_2 ]
    
    print "plot..."
    plot_multiway_tree("MultiwayFunctionTree[%d-%d][%d-%d][%s>%f].png"%(range1+range2+("gcai",gcai_threshold)),
                       DF_EDGE_LIST, FF_EDGE_LIST, FA_EDGE_LIST,
                       DOMAIN_SET_1, FUNC_SET_1, ANCE_SET_1,
                       DOMAIN_SET_2, FUNC_SET_2, ANCE_SET_2,
                       func_commun, ance_commun)

def save_key_vars(fname, func_meta_list, func_op_list,
                  DF_EDGE_LIST_meta, FF_EDGE_LIST_meta, FA_EDGE_LIST_meta,
                  DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta,                  
                  DF_EDGE_LIST_op, FF_EDGE_LIST_op, FA_EDGE_LIST_op,
                  DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op,                  
                  DF_EDGE, FF_EDGE, FA_EDGE,
                  base_commun, func_commun, ance_commun ):

    handler = open(fname,'w')
    
    tvDict  = {"func_meta_list":func_meta_list,
               "func_op_list":func_op_list,
               "DF_EDGE_LIST_meta":DF_EDGE_LIST_meta,
               "FF_EDGE_LIST_meta":FF_EDGE_LIST_meta,
               "FA_EDGE_LIST_meta":FA_EDGE_LIST_meta,
               "DOMAIN_SET_meta":DOMAIN_SET_meta,
               "FUNC_SET_meta":FUNC_SET_meta,
               "ANCE_SET_meta":ANCE_SET_meta,
               "DF_EDGE_LIST_op":DF_EDGE_LIST_op,
               "FF_EDGE_LIST_op":FF_EDGE_LIST_op,
               "FA_EDGE_LIST_op":FA_EDGE_LIST_op,
               "DOMAIN_SET_op":DOMAIN_SET_op,
               "FUNC_SET_op":FUNC_SET_op,
               "ANCE_SET_op":ANCE_SET_op,
               "DF_EDGE":DF_EDGE,
               "FF_EDGE":FF_EDGE,
               "FA_EDGE":FA_EDGE,
               "base_commun":base_commun,
               "func_commun":func_commun,
               "ance_commun":ance_commun}

    cPickle.dump(tvDict, handler, 2)
    handler.close()
    
def load_key_vars(fname):
    handler = open(fname)
    tvDict = cPickle.load(handler)
    handler.close()
    return ( tvDict["func_meta_list"], tvDict["func_op_list"],
    tvDict["DF_EDGE_LIST_meta"], tvDict["FF_EDGE_LIST_meta"], tvDict["FA_EDGE_LIST_meta"],
    tvDict["DOMAIN_SET_meta"], tvDict["FUNC_SET_meta"], tvDict["ANCE_SET_meta"],
    tvDict["DF_EDGE_LIST_op"], tvDict["FF_EDGE_LIST_op"], tvDict["FA_EDGE_LIST_op"],
    tvDict["DOMAIN_SET_op"], tvDict["FUNC_SET_op"], tvDict["ANCE_SET_op"],
    tvDict["DF_EDGE"], tvDict["FF_EDGE"], tvDict["FA_EDGE"],
    tvDict["base_commun"], tvDict["func_commun"], tvDict["ance_commun"])

def save_edges_to_csv(fname, func_edge_list, label_dict) :
    handler = open(fname,'w')
    writer = csv.writer(handler)
    for sf,tf in func_edge_list :
        writer.writerow([sf,tf,'pp','TRUE',label_dict[sf], label_dict[tf] ])
    handler.close()

    
def plot_sub_2(rangeCeiling=0,rangeFloor=50,gcai_threshold=0.6,tvfname = "/tmp/.tmpsave"):
    
    print "loading..."
    # pfam2go_dict, goslimmeta_dict, goDict, cais_dict_meta, domain_abundance_meta = load_all_metagenomic()
    # pfam2go_dict, goslimmeta_dict, goDict, cais_dict_meta, domain_abundance_meta = load_all_meta_by_meta()
    pfam2go_dict, goslimmeta_dict, goDict, cais_dict_meta, domain_abundance_meta = load_all_meta_by_meta_moy()
    cais_dict_op, domain_abundance_op = load_all_opposite()
    
    print "convert to fonction domain list..."
    go_ab_meta_dict = fonction_abondance(pfam2go_dict, domain_abundance_meta)
    go_gcai_meta_dict = fonction_max_gcai(pfam2go_dict, cais_dict_meta)
    go_ab_op_dict = fonction_abondance(pfam2go_dict, domain_abundance_op)
    go_gcai_op_dict = fonction_max_gcai(pfam2go_dict, cais_dict_op)

    print "filtering by gcai threshold..."
    go_ab_meta_dict_clean = filter_go_ab_dict_by_gcai(go_ab_meta_dict, go_gcai_meta_dict, gcai_threshold)
    go_ab_op_dict_clean = filter_go_ab_dict_by_gcai(go_ab_op_dict, go_gcai_op_dict, gcai_threshold)

    print "[Statistic] Meta group has %d sequences over %f gcai\n"%(len(go_ab_meta_dict_clean.keys()),gcai_threshold)
    if len(go_ab_meta_dict_clean.keys()) < rangeFloor :
        print "Number of metagenomic funtion is less than max of range requist."
        exit()
        
    print "[Statistic] Oppo group has %d sequences over %f gcai\n"%(len(go_ab_op_dict_clean.keys()),gcai_threshold)
    if len(go_ab_op_dict_clean.keys()) < rangeFloor :
        print "Number of opposite funtion is less than max of range requist."
        exit()
    
    print "sort functions by abondance..."
    go_ab_meta_list = sortDictByValue(go_ab_meta_dict_clean)
    func_meta_list  = np.array(go_ab_meta_list)[:,0][rangeCeiling:rangeFloor]
    go_ab_op_list = sortDictByValue(go_ab_op_dict_clean)
    func_op_list  = np.array(go_ab_op_list)[:,0][rangeCeiling:rangeFloor]

    print "find all parents of metagenome"
    reset_global()
    find_all_fonction_edges(func_meta_list, goslimmeta_dict, goDict)
    DF_EDGE_LIST_meta, FF_EDGE_LIST_meta, FA_EDGE_LIST_meta = copy_global_list()
    DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta = copy_global_set()
    
    print "find all parents of opposite"
    reset_global()
    find_all_fonction_edges(func_op_list, goslimmeta_dict, goDict)
    DF_EDGE_LIST_op, FF_EDGE_LIST_op, FA_EDGE_LIST_op = copy_global_list()    
    DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op = copy_global_set()

    print "Merge all edges..."
    def merge_list(list1,list2):
        lo = list1
        for l in list2:
            if not l in list1 :
                lo.append(l)
        return lo
        
    DF_EDGE = merge_list(DF_EDGE_LIST_meta,DF_EDGE_LIST_op)
    FF_EDGE = merge_list(FF_EDGE_LIST_meta,FF_EDGE_LIST_op)
    FA_EDGE = merge_list(FA_EDGE_LIST_meta,FA_EDGE_LIST_op)
    
    print "searching commun functions between metagenome et opposite ..."
    base_commun = [ f for f in func_meta_list if f in func_op_list ]
    func_commun = [ f for f in FUNC_SET_meta if f in FUNC_SET_op ]
    ance_commun = [ f for f in ANCE_SET_meta if f in ANCE_SET_op ]

    print "saving variances..."
    save_key_vars(tvfname, func_meta_list, func_op_list,
                  DF_EDGE_LIST_meta, FF_EDGE_LIST_meta, FA_EDGE_LIST_meta,
                  DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta,                  
                  DF_EDGE_LIST_op, FF_EDGE_LIST_op, FA_EDGE_LIST_op,
                  DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op,                  
                  DF_EDGE, FF_EDGE, FA_EDGE,
                  base_commun, func_commun, ance_commun )

    '''
    print "plot..."
    plot_multiway_tree("MultiwayFunctionTreeAllMMoyOppo[%d-%d][%s>%f].png"%(rangeCeiling,rangeFloor,"gcai",gcai_threshold),
                       DF_EDGE, FF_EDGE, FA_EDGE,
                       DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta,
                       DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op,
                       func_commun, ance_commun,base_commun,
                       func_meta_list, func_op_list )
    '''
    
def plot_fast_sub_2(rangeCeiling=0,rangeFloor=50,gcai_threshold=0.6,tvfname = "/tmp/.tmpsave"):
    print "loading variances..."
    ( func_meta_list, func_op_list,
    DF_EDGE_LIST_meta, FF_EDGE_LIST_meta, FA_EDGE_LIST_meta,
    DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta,                  
    DF_EDGE_LIST_op, FF_EDGE_LIST_op, FA_EDGE_LIST_op,
    DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op,                  
    DF_EDGE, FF_EDGE, FA_EDGE,
    base_commun, func_commun, ance_commun ) = load_key_vars(tvfname)

    goDict = load_godict()
    label_dict = {}
    label_func_list = func_meta_list.tolist() + func_op_list.tolist() + func_commun + ance_commun + base_commun + list(ANCE_SET_meta) + list(ANCE_SET_op) + list(FUNC_SET_meta) + list(FUNC_SET_op)
    
    for f in label_func_list :
        label_dict[f] = goDict[f].name
    '''    
        label_dict[f] = f + "\n" + goDict[f].name
    print "plot..."
    plot_multiway_tree("MultiwayFunctionTreeAllMMoyOppo[%d-%d][%s>%.2f].png"%(rangeCeiling,rangeFloor,"gcai",gcai_threshold),
                       DF_EDGE, FF_EDGE, FA_EDGE,
                       DOMAIN_SET_meta, FUNC_SET_meta, ANCE_SET_meta,
                       DOMAIN_SET_op, FUNC_SET_op, ANCE_SET_op,
                       func_commun, ance_commun,base_commun,
                       func_meta_list, func_op_list,
                       label_func_list, label_dict)
    '''
    print "saving edges in csv format..."
    fname = "MultiwayFunctionTreeAllMMoyOppoEdges[%d-%d][%s>%.2f].csv"%(rangeCeiling,rangeFloor,"gcai",gcai_threshold)
    save_edges_to_csv(fname, DF_EDGE_LIST_meta + FF_EDGE_LIST_meta + FA_EDGE_LIST_meta + DF_EDGE_LIST_op + FF_EDGE_LIST_op + FA_EDGE_LIST_op, label_dict)
    
# plot_sub_1((0,5),(5,10))
# plot_sub_1((0,10),(10,20))
# plot_sub_1((0,50),(50,100))
# plot_sub_2(rangeCeiling=0,rangeFloor=50)
plot_fast_sub_2(rangeCeiling=0,rangeFloor=50)

