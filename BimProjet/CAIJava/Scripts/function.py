#!env python

import pydot
import cPickle
import numpy as np
# import networkx as nx
from graphviz import Digraph

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
        
