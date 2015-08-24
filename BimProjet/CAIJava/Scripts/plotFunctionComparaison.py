#!env python
import csv
import numpy as np
import operator

from function import loadToDict, load_pfam2go_toDict, GOAnnotation, plot_switch_view_compare, switch_view, domain_function_list, domain_function_dict, familyReference
from BimProjetLib import readClstr, readgCAIs, read2step, getGIDList, domainIdDict, domainGenomeDict, accumulateNivExprDom, sortDictByValue, gCAIsDictToDomDict

with open("pfam2go") as handle :
    pfam2go_dict = load_pfam2go_toDict(handle)

with  open("GO_SLIM_META_DICT.txt",'r') as handle :
    goslimmeta_dict = loadToDict(handle)

with open("GO_DICT.txt",'r') as handle :
    goDict = loadToDict(handle)

############################################################
def sub_1(vseuil=0.9):
    # loading
    CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
    # CAI_dict = readgCAIs("../output/cais.lst.2step")
    CAI_dict = readgCAIs("../output/cais_len30.lst")
    step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
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
    deSortList,dg,de = familyReference(nivExprAccuDomDict,gCAIsDomDict,domains,pfam2go_dict,neSeuil=0,gSeuil=0)
    klist = []
    ndg = {}
    nde = {}
    for k,vl in dg.iteritems():
        if np.max(vl) > vseuil:
            klist.append(k)
    for k in klist:
        ndg[k] = dg[k]
        nde[k] = de[k]
    deSortList = np.array(sortDictByValue(nde))
    dfList = domain_function_list(np.array(deSortList)[:,0],pfam2go_dict)
    bio_process_1, molec_func_1, cellu_comp_1 = switch_view(dfList,goDict,goslimmeta_dict,ndg,nde)
    return bio_process_1, molec_func_1, cellu_comp_1
    
############################################################
def sub_2():
    dom_gcai = {}
    dom_abundance = {}
    with open("/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/CAIJava_pipeline/domain_gcai_abundance.csv") as handle :
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
    da = sorted(dom_abundance.items(), key=operator.itemgetter(1),reverse=True)
    dom_list = np.array(da)[:,0]
    dfList = domain_function_list(dom_list,pfam2go_dict)
    bio_process_2,molec_func_2,cellu_comp_2 = switch_view(dfList,goDict,goslimmeta_dict,dom_gcai,dom_abundance)
    return bio_process_2,molec_func_2,cellu_comp_2
##############################################################

bio_process_1, molec_func_1, cellu_comp_1 = sub_1(0.6)
bio_process_2, molec_func_2, cellu_comp_2 = sub_2()

# plot_switch_view_compare(func_type,func_dict_1,func_dict_2,fname=None,figsize=None,caiSeuil=0.0)
plot_switch_view_compare('Biological Process',bio_process_1,bio_process_2,'BiologicalProcessCmp.png',(70,80))
plot_switch_view_compare('Cellular Component',cellu_comp_1,cellu_comp_2, 'CellularComponentCmp.png',(70,30))
plot_switch_view_compare('Molecular Function',molec_func_1, molec_func_2, 'MolecularFunctionCmp.png',(70,90))
