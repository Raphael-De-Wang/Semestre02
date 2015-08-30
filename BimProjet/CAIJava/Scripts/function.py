#!env python

import csv
import pydot
import cPickle
import operator
import numpy as np
# import networkx as nx
from graphviz import Digraph
import matplotlib.pyplot as plt

from BimProjetLib import sortDictByValue

def loadToDict(handle):
    goDict = cPickle.load(handle)
    return goDict
    
class GOAnnotation():
    def __init__(self,id,acc,name,type,parent=None):
        self.id               = id;
        self.count            = 0;
        self.countNorm        = 0;
        self.totalCount       = 0;
        self.totalCountNorm   = 0;
        self.acc              = acc;
        self.name             = name;
        self.type             = type;
        self.clusterDist      = [];
        self.domains          = {};
        self.sequences        = {};
        self.descendants      = [];
        self.parent           = [];
        self.firstDescendants = [];
    
    def addDescendant(self,acc):
        print type(acc)
        if not acc == self.acc and not acc in self.descendants:
            self.descendants.append(acc);
        
    def addFirstDescendant(self,acc):
        if not acc == self.acc and not acc in self.firstDescendants:
            self.firstDescendants.append(acc);

def load_pfam2go_toDict(handle):
    pfam2go = dict()
    for line in handle:
        if line[:4] == 'Pfam':
            title, context = line.split(' > ')
            domain,tag = title.split()
            desc, func = context.split(' ; ')
            func = func.strip()
            domain = domain[5:]
            if not pfam2go.has_key(domain):
                pfam2go[domain] = [[tag,desc,func]]
            else:
                pfam2go[domain].append([tag,desc,func])
    return pfam2go
    
def pfam2go_to_funcDomainDict(pfam2go):
    fd_dict = dict()
    for domain in pfam2go:
        for tag,desc,func in pfam2go[domain]:
            if not fd_dict.has_key(func):
                fd_dict[func] = []
            if domain not in fd_dict[func]:
                fd_dict[func].append(domain)
    return fd_dict

DOMAIN_LIST = []
FUNC_LIST   = []
DF_LIST     = []
FF_LIST     = []

def add_df_edge(dot,domain,funcName):
    func = funcName.replace(':','')
    if [domain,funcName] not in DF_LIST:
        if domain not in DOMAIN_LIST:
            # node1 = pydot.Node(domain, style="filled", fillcolor="green")
            node1 = pydot.Node(domain, style="filled", fillcolor="green",layer='1',nodesep="1.5")
            DOMAIN_LIST.append(domain)
            dot.add_node(node1)
        else:
            node1 = dot.get_node(domain)
            if type(node1) == list:
                node1 = node1[0]
        #####
        if funcName not in FUNC_LIST:
            node2 = pydot.Node(func, style="filled", fillcolor="yellow", layer='2',nodesep="1.5")
            FUNC_LIST.append(funcName)
            dot.add_node(node2)
        else:
            node2 = dot.get_node(func)
            if type(node2) == list:
                node2 = node2[0]
        #####
        edge = pydot.Edge(node1,node2)
        DF_LIST.append([domain,func])
        dot.add_edge(edge)

def add_da_edge(dot,domain,anceName):
    ance = anceName.replace(':','')
    if [domain,anceName] not in DF_LIST:
        if domain not in DOMAIN_LIST:
            # node1 = pydot.Node(domain, style="filled", fillcolor="green")
            node1 = pydot.Node(domain, style="filled", fillcolor="green",layer='1',nodesep="1.5")
            DOMAIN_LIST.append(domain)
            dot.add_node(node1)
        else:
            node1 = dot.get_node(domain)
            if type(node1) == list:
                node1 = node1[0]
        #####
        if anceName not in FUNC_LIST:
            node2 = pydot.Node(ance, style="filled", fillcolor="#976856", layer='3',nodesep="1.5")
            FUNC_LIST.append(anceName)
            dot.add_node(node2)
        else:
            node2 = dot.get_node(ance)
            if type(node2) == list:
                node2 = node2[0]
        #####
        edge = pydot.Edge(node1,node2)
        DF_LIST.append([domain,anceName])
        dot.add_edge(edge)
        
def add_ff_edge( dot, funcName1, funcName2 ):
    func1 = funcName1.replace(':','')
    func2 = funcName2.replace(':','')
    if [funcName1,funcName2] not in FF_LIST:
        if funcName1 not in FUNC_LIST:
            node1 = pydot.Node(func1, style="filled", fillcolor="yellow", layer='2',nodesep="1.5")
            FUNC_LIST.append(funcName1)
            dot.add_node(node1)
        else:
            node1 = dot.get_node(func1)
            if type(node1) == list:
                node1 = node1[0]
        #####
        if funcName2 not in FUNC_LIST:
            node2 = pydot.Node(func2, style="filled", fillcolor="yellow", layer='2',nodesep="1.5")
            FUNC_LIST.append(funcName2)
            dot.add_node(node2)
        else:
            node2 = dot.get_node(func2)
            if type(node2) == list:
                node2 = node2[0]
        #####
        edge = pydot.Edge(node1,node2)
        FF_LIST.append([funcName1,funcName2])
        dot.add_edge(edge)

def add_fa_edge( dot, funcName, anceName ):
    func = funcName.replace(':','')
    ance = anceName.replace(':','')
    if [funcName,anceName] not in FF_LIST:
        if funcName not in FUNC_LIST:
            node1 = pydot.Node(func, style="filled", fillcolor="yellow", layer='1',nodesep="1.5")
            FUNC_LIST.append(funcName)
            dot.add_node(node1)
        else:
            node1 = dot.get_node(func)
            if type(node1) == list:
                node1 = node1[0]
        #####
        if anceName not in FUNC_LIST:
            node2 = pydot.Node(ance, style="filled", fillcolor="#976856", layer='3',nodesep="1.5")
            FUNC_LIST.append(anceName)
            dot.add_node(node2)
        else:
            node2 = dot.get_node(ance)
            if type(node2) == list:
                node2 = node2[0]
        #####
        edge = pydot.Edge(node1,node2)
        FF_LIST.append([funcName,anceName])
        dot.add_edge(edge)
        
def search_ancestor( dot, orig, funcName, go_dict, meta_dict, count=0):
    if count > 10:
        # print count, 'times call search ancestor'
        return count
    for key in go_dict:
        if funcName in go_dict[key].descendants:
            # print orig,'->',funcName
            if key in meta_dict.keys():
                add_fa_edge( dot, funcName, key )
            else:
                add_ff_edge( dot, funcName, key )
                count += 1
                count = search_ancestor( dot, funcName, key, go_dict, meta_dict, count)
    return count
        
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

def plot_switch_view(func_type,func_dict,fname=None,figsize=None,caiSeuil=0.6):
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
        cais_list = [ j for i in cais_list for j in i ]        
        if np.max(cais_list) < caiSeuil :
            continue
        if not fndict.has_key(key):
            fndict[key] = sum(ne_list)
        else :
            fndict[key] += sum(ne_list)
    for key,mNe in sortDictByValue(fndict,False):
        value = func_dict[key]
        (name,dom_list,cais_list,ne_list) = value
        cais_list = [ j for i in cais_list for j in i ]
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
            width=domList, bottom=bottoms,align="edge", # "center",
            color="#2aa198",orientation="horizontal",height=1.0)
    plt.yticks(bottoms+0.1,nameList)
    ax = plt.subplot(132,axisbg="#fdf6e3")
    ax.set_title("Niveau d'Expression")
    plt.bar(left=np.zeros(len(nameList)),
            width=neList, bottom=bottoms,align="edge", # "center",
            color="#2aa198",orientation="horizontal",height=1.0)
    plt.yticks(bottoms+0.1,funcList)
    ax = plt.subplot(133,axisbg="#fdf6e3")
    ax.set_title("gCAI")
    plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMaxList, bottom=bottoms,
            color="b",orientation="horizontal",
            height=0.7,label='Max')
    plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMeanList, bottom=bottoms,
            color="g",orientation="horizontal",
            height=0.8,label='Mean')
    plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMinList, bottom=bottoms,
            color="r",orientation="horizontal",
            height=0.9,label='Min')
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    plt.yticks(bottoms+0.1,funcList)
    # ax.yaxis.tick_right()
    # fig.subplots_adjust(left=0.3,hspace=0.05,bottom=0.1,top=0.90)
    fig.subplots_adjust(left=0.3,bottom=0.1,top=0.90)
    if fname == None:
        plt.show()
    else:
        plt.savefig(fname)
        '''
        handle = open(fname.replace('png','txt'),'w')
        for f in funcList:
            handle.write("%s\n"%f)
        handle.close()
        '''
    plt.close(fig)

def plot_switch_view_all_in_one(bio_process_1, bio_process_2, molec_func_1,molec_func_2 , cellu_comp_1, cellu_comp_2, fname=None, figsize=None, caiSeuil=0.0):
    def plot_switch_view(func_type,func_dict_1,func_dict_2,fname,caiSeuil,subplot_pos):
        nameList = []
        domList  = []
        caisMinList = []
        caisMaxList = []
        caisMeanList= []
        caisMinList_2 = []
        caisMaxList_2 = []
        caisMeanList_2= []
        neList      = []
        funcList    = []
        fndict      = {}
        for key,value in func_dict_1.iteritems():
            (name,dom_list,cais_list,ne_list) = value
            (name_2,dom_list_2,cais_list_2,ne_list_2) = func_dict_2[key]
            cais_list = [ j for i in cais_list for j in i ]        
            cais_list_2 = [ j for i in cais_list_2 for j in i ]
            if name <> name_2 or len(cais_list) <> len(cais_list_2) or len(dom_list) <> len(dom_list_2) or len(ne_list) <> len(ne_list_2) :
                raise ValueError("Data Conflict In func_dict")
            if np.max(cais_list) < caiSeuil :
                continue
            if not fndict.has_key(key):
                fndict[key] = sum(ne_list)
            else :
                fndict[key] += sum(ne_list)
        for key,mNe in sortDictByValue(fndict,False):
            value_1 = func_dict_1[key]
            value_2 = func_dict_2[key]
            (name,dom_list,cais_list,ne_list) = value_1
            (name_2,dom_list_2,cais_list_2,ne_list_2) = value_2
            cais_list = [ j for i in cais_list for j in i ]
            cais_list_2 = [ j for i in cais_list_2 for j in i ]
            caisMinList.append(min(cais_list))
            caisMaxList.append(max(cais_list))
            caisMeanList.append(np.mean(cais_list))
            caisMinList_2.append(min(cais_list_2))
            caisMaxList_2.append(max(cais_list_2))
            caisMeanList_2.append(np.mean(cais_list_2))
            nameList.append(name)
            domList.append(len(dom_list))
            neList.append(sum(ne_list))
            funcList.append(key)
        bottoms = np.arange(len(nameList))*1.2
        ax1 = plt.subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2]+1,axisbg="#fdf6e3")
        ax1.set_title("Nombre de Domain par Function")
        plt.bar(left=np.zeros(len(nameList)),
            width=domList, bottom=bottoms,align="edge", # "center",
            color="#2aa198",orientation="horizontal",height=1.0)
        plt.yticks(bottoms+0.1,nameList)
        ax2 = plt.subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2]+2,axisbg="#fdf6e3")
        ax2.set_title("Niveau d'Expression")
        plt.bar(left=np.zeros(len(nameList)),
            width=neList, bottom=bottoms,align="edge", # "center",
            color="#2aa198",orientation="horizontal",height=1.0)
        plt.yticks(bottoms+0.1,funcList)
        ax3 = plt.subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2]+3,axisbg="#fdf6e3")
        ax3.set_title("gCAI")
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMaxList, bottom=bottoms,
            color="b",orientation="horizontal",
            height=0.7)
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMeanList, bottom=bottoms,
            color="g",orientation="horizontal",
            height=0.8)
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMinList, bottom=bottoms,
            color="r",orientation="horizontal",
            height=0.9)
        plt.yticks(bottoms+0.1,funcList)
        ax4 = plt.subplot(subplot_pos[0],subplot_pos[1],subplot_pos[2]+4,axisbg="#fdf6e3")
        ax4.set_title("gCAI 2")
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMaxList_2, bottom=bottoms,
            color="b",orientation="horizontal",
            height=0.7,label='Max')
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMeanList_2, bottom=bottoms,
            color="g",orientation="horizontal",
            height=0.8,label='Mean')
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
            width=caisMinList_2, bottom=bottoms,
            color="r",orientation="horizontal",
            height=0.9,label='Min')
        plt.yticks(bottoms+0.1,funcList)
        return (ax1,ax2,ax3,ax4)
    fig = plt.figure(figsize=figsize)
    # plt.suptitle(, fontsize=15)
    ax1,ax2,ax3,ax4 = plot_switch_view("Biological Process",bio_process_1,bio_process_2,fname,caiSeuil,[3,4,0])
    bb = ax1.get_position()
    plt.figtext(.1,bb.max[1]-.05,'Biological Process', fontsize=40, ha='center')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    ax1,ax2,ax3,ax4 = plot_switch_view("Molecular Function",molec_func_1,molec_func_2,fname,caiSeuil, [3,4,4])
    bb = ax1.get_position()
    plt.figtext(.1,bb.max[1]-.05,'Molecular Function',fontsize=40, ha='center')
    ax1,ax2,ax3,ax4 = plot_switch_view("Cellular Component",cellu_comp_1,cellu_comp_2,fname,caiSeuil, [3,4,8])
    bb = ax1.get_position()
    plt.figtext(.1,bb.max[1]-0.05,'Cellular Component',fontsize=40, ha='center')
    fig.subplots_adjust(left=0.3,bottom=0.1,top=0.90)
    if fname == None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)

def plot_switch_view_compare(func_type,func_dict_1,func_dict_2,fname=None,figsize=None,caiSeuil=0.0):
    nameList = []
    funcList = []
    key_list = []    
    fndict_1 = {}
    for key,value_1 in func_dict_1.iteritems():
        (name_1,dom_list_1,cais_list_1,ne_list_1) = value_1
        if max(cais_list_1) < caiSeuil :
            continue
        if func_dict_2.has_key(key):
            if not fndict_1.has_key(key):
                fndict_1[key] = sum(ne_list_1)
            else :
                fndict_1[key] += sum(ne_list_1)
    for key,mNe in sortDictByValue(fndict_1,False):
        key_list.append(key)
        (name,dom_list,cais_list,ne_list) = func_dict_1[key]
        nameList.append(name)
        funcList.append(key)
    bottoms = np.arange(len(nameList))*1.2
    fig = plt.figure(figsize=figsize)
    plt.suptitle(func_type, fontsize=15)
    def sub_plot_switch_view(func_dict, key_list, sub_num, name_appendix):
        domList  = []
        neList   = []
        caisMinList = []
        caisMaxList = []
        caisMeanList= []
        for key in key_list:
            value = func_dict[key]
            (name,dom_list,cais_list,ne_list) = value
            cais_list = [j for i in cais_list for j in i ]
            caisMinList.append(min(cais_list))
            caisMaxList.append(max(cais_list))
            caisMeanList.append(np.mean(cais_list))
            domList.append(len(dom_list))
            neList.append(sum(ne_list))
        if (sub_num==1):
            ax = plt.subplot(161,axisbg="#fdf6e3")
        else:
            ax = plt.subplot(162,axisbg="#fdf6e3")
        ax.set_title("Nombre de Domain par Function %s"%name_appendix)
        plt.bar(left=np.zeros(len(nameList)),
                width=domList, bottom=bottoms,align="edge", # "center",
                color="#2aa198",orientation="horizontal",height=1.0)
        plt.yticks(bottoms+0.1,nameList)
        if (sub_num==1):
            ax = plt.subplot(163,axisbg="#fdf6e3")
        else:
            ax = plt.subplot(164,axisbg="#fdf6e3")
        ax.set_title("Niveau d'Expression %s"%name_appendix)
        plt.bar(left=np.zeros(len(nameList)),
                width=neList, bottom=bottoms,align="edge", # "center",
                color="#2aa198",orientation="horizontal",height=1.0)
        plt.yticks(bottoms+0.1,funcList)
        if (sub_num==1):        
            ax = plt.subplot(165,axisbg="#fdf6e3")
        else:
            ax = plt.subplot(166,axisbg="#fdf6e3")
            
        ax.set_title("gCAI %s"%name_appendix)
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
                width=caisMaxList, bottom=bottoms,
                color="b",orientation="horizontal",
                height=0.7,label='Max')
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
                width=caisMeanList, bottom=bottoms,
                color="g",orientation="horizontal",
                height=0.8,label='Mean')
        plt.bar(left=np.zeros(len(nameList)),align="edge", # "center",
                width=caisMinList, bottom=bottoms,
                color="r",orientation="horizontal",
                height=0.9,label='Min')
        # plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        plt.yticks(bottoms+0.1,funcList)
        # ax.yaxis.tick_right()
        # fig.subplots_adjust(left=0.3,hspace=0.05,bottom=0.1,top=0.90)
    sub_plot_switch_view(func_dict_1, key_list, 1, "Metagenomic")
    sub_plot_switch_view(func_dict_2, key_list, 2, "")
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

def load_pfam2go():
    with open("pfam2go") as handle :
        pfam2go_dict = load_pfam2go_toDict(handle)
    return pfam2go_dict
    
def load_goslim():        
    with  open("GO_SLIM_META_DICT.txt",'r') as handle :
        goslimmeta_dict = loadToDict(handle)
    return goslimmeta_dict
        
def load_godict():        
    with open("GO_DICT.txt",'r') as handle :
        goDict = loadToDict(handle)
    return goDict

def load_domain_gcai_abundance(fname="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline/domain_gcai_abundance.csv"):
    dom_gcai = {}
    dom_abundance = {}
    with open(fname) as handle :
        reader = csv.reader(handle, delimiter='\t')
        for dname, gcai, abundance in reader :
            if not dom_abundance.has_key(dname) :
                dom_abundance[dname] = int(abundance)
            elif dom_abundance[dname] <> int(abundance) :
                raise ValueError("[%s] Abundance Value Conflict. "%dname)
            if not dom_gcai.has_key(dname) :
                dom_gcai[dname] = [float(gcai)]
            else:
                dom_gcai[dname].append(float(gcai))
    return (dom_gcai, dom_abundance)

def load_all_func_data():
    pfam2go_dict = load_pfam2go()
    goslimmeta_dict = load_goslim()
    goDict = load_godict()
    (dom_gcai, dom_abundance) = load_domain_gcai_abundance()
    return ( pfam2go_dict, goslimmeta_dict, goDict, dom_gcai, dom_abundance )

