#!env python

import csv
import sys
import pylab
import operator
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from BimProjetLib import readClstr, readgCAIs, read2step, getGIDList, domainIdDict, domainGenomeDict, accumulateNivExprDom, sortDictByValue, gCAIsDictToDomDict

def plotgCAIsDomain(X,Y,xDiv=1.,yDiv=20.,figName=None):
    fig = pylab.figure(figsize=(320,3))
    # pylab.hist2d(X,Y,bins=[len(X)/xDiv,yDiv], norm=mpl.colors.LogNorm(), cmap=mpl.cm.jet)
    pylab.hist2d(X,Y,bins=[len(np.unique(X))/xDiv,yDiv],norm=mpl.colors.LogNorm(), cmap=mpl.cm.jet)
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.title("Desity Plot - Domains Grouped by [%d], gCAIs Divided by [%d], Top [%d] Domains (PS : values in brackets are manipulable)"%(xDiv,yDiv,len(np.unique(Y))))
    fig.subplots_adjust(bottom=0.25)
    cbar = pylab.colorbar()
    if figName:
        pylab.savefig(figName)

    else:
        pylab.show()
    pylab.close(fig)

def plotgCAIsDomain2(X,Y1,Y2,figName=None,xticks=[]):
    T = np.arctan(Y2) #,np.ones(len(X)))
    fig = pylab.figure(figsize=(30,8))
    xUniLen = len(np.unique(X))
    pylab.scatter(X, Y1, s=5, c=T, alpha=.4, marker='x')
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.ylim(0,1.02), pylab.yticks(pylab.arange(0,1.1,0.1))
    pylab.xlim(0,xUniLen), pylab.xticks(xticks)
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen)*0.025,'k-')
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen),'k-')
    pylab.title("gCAIs values distributions ordered by expression level, the colors correspond to gCAI value, each colume represent a domain")
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close()
    
def plotgCAIsDomainCompare(X,Y1,Y2,figName=None,xticks=[]):
    T1 = np.arctan(Y1) #,np.ones(len(X)))
    T2 = np.arctan(Y2) #,np.ones(len(X)))
    fig = pylab.figure(figsize=(60,8))
    xUniLen = len(np.unique(X))
    pylab.scatter(X, Y1, s=5, c=T1, alpha=.4, marker='x')
    pylab.scatter(X, Y1*(-1.), s=5, c=T2, alpha=.4, marker='x')
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.ylim(-1.02,1.02),pylab.yticks(np.array(range(10,-11,-1))/10.)
    # pylab.yticks(pylab.arange(0,1.1,0.1).tolist()+pylab.arange(1.,0,-0.1).tolist())
    pylab.xlim(0,xUniLen), pylab.xticks(xticks)
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen)*(-1.),'k-')
    pylab.plot(np.arange(xUniLen),np.zeros(xUniLen),'k-')
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen),'k-')
    pylab.plot([50 for i in range(21)],np.arange(-1,1.1,0.1),'k-')
    pylab.plot([100 for i in range(21)],np.arange(-1,1.1,0.1),'k-')
    pylab.plot([200 for i in range(21)],np.arange(-1,1.1,0.1),'k-')
    pylab.plot([300 for i in range(21)],np.arange(-1,1.1,0.1),'k-')
    pylab.plot([400 for i in range(21)],np.arange(-1,1.1,0.1),'k-')
    pylab.title("gCAIs values distributions ordered by expression level, the colors correspond to gCAI value, each colume represent a domain")
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close()
    
def load_ORF_gcai(handle):
    ORF_dict = {}
    for line in handle:
        line = line.split()
        if line[1] == 'CAI':
            continue
        ORF_dict[line[0]] = float(line[1])
    return ORF_dict
            
def sub_load_sorted_by_abundance():
    # loading
    CLS_dict = readClstr("AT_arc_metatrans.filtered.fasta.clstr")
    step2MList = read2step("AT_arc_metatrans.filtered.fasta.6RF.faa.e_minus10_pos_neg_covSeq_EM10.archs.2step")
    gIDList = getGIDList(step2MList)
    # Accumulate Niveau Expression :
    domGenomeDict = domainGenomeDict(gIDList)
    nivExprAccuDomDict = accumulateNivExprDom(CLS_dict,domGenomeDict)
    nivExprAccuDomSortList = np.array(sortDictByValue(nivExprAccuDomDict))
    return nivExprAccuDomSortList

def build_domain_orf_dict(orf_list):
    domain_orf_dict = {}
    for orf in orf_list:
        domain = orf.split('_')[-1]
        if domain_orf_dict.has_key(domain):
            domain_orf_dict[domain].append(orf)
        else:
            domain_orf_dict[domain]=[orf]
    return domain_orf_dict
    
def build_domain_orf1_orf2_list(domain_list, domain_orf_dict, orf_dict_1, orf_dict_2):
    dID = 0
    dlist = []
    orf_list = []
    dID_orf1_orf2_list = []
    for d in domain_list :
        dID += 1
        if not domain_orf_dict.has_key(d):
            continue
        for orf in domain_orf_dict[d]:
            try:
                dID_orf1_orf2_list.append([dID,orf_dict_1[orf],orf_dict_2[orf]])
                dlist.append(d)
                orf_list.append(orf)
            except:
                print orf_dict_1[orf]
                # raise
    return np.array(dID_orf1_orf2_list),dlist,orf_list
    
def persist_domain_orf1_orf2_list(label_list, orf1_list, orf2_list, outpath):
    with open(outpath,'w') as output_handle:
        for ind,label in enumerate(label_list) :
            output_handle.write('\t'.join([label,'%f'%orf1_list[ind],'%f'%orf2_list[ind]])+'\n')

def load_domain_orf1_orf2_list(inpath):
    label_list = []
    orf1_list  = []
    orf2_list  = []
    with open(inpath) as input_handle:
        for line in input_handle:
            label,orf1,orf2 = line.strip().split('\t')
            label_list.append(label)
            orf1_list.append(float(orf1))
            orf2_list.append(float(orf2))
    return (label_list, orf1_list, orf2_list)

def sub_1(inum):
    ifname_1 = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_len%d.lst"%inum
    ifname_2 = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/cais_ewvalue_len%d.lst"%inum
    ofname   = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/domain_cais_1_2_len%d.lst"%inum
    handle_1 = open(ifname_1)
    handle_2 = open(ifname_2)
    orf_dict_1 = load_ORF_gcai(handle_1)
    orf_dict_2 = load_ORF_gcai(handle_2)
    domain_orf_dict = build_domain_orf_dict(orf_dict_1.keys())
    nivExprAccuDomSortList = sub_load_sorted_by_abundance()
    d12,dlist,orf_list = build_domain_orf1_orf2_list(nivExprAccuDomSortList[:,0], domain_orf_dict, orf_dict_1, orf_dict_2)
    orf_list, d12 = sort_domain_orf1_orf2_list(orf_list, d12.tolist())
    plotgCAIsDomainCompare(d12[:,0],d12[:,1],d12[:,2],figName="gCAIsDistrColorBy2_Compare_len%d"%inum)
    # persist_domain_orf1_orf2_list(orf_list, d12[:,1], d12[:,2], ofname)
    handle_1.close()
    handle_2.close()

def filter_gcais_less_important(label_list, orf1_list, orf2_list, seuil=0.9):
    domain_orf1_orf2_list = []
    for ind, label in enumerate(label_list) :
        if orf1_list[ind] >= seuil or orf2_list[ind] >= seuil :
            domain_orf1_orf2_list.append([label, orf1_list[ind], orf2_list[ind], abs(orf1_list[ind] - orf2_list[ind])])
    return domain_orf1_orf2_list

def sort_domain_orf1_orf2_list(orf_list, d12_list):
    all_in_list = [ [orf_list[i]] + d12_list[i] for i in range(len(orf_list)) ]
    all_in_list = sorted(all_in_list, key=operator.itemgetter(3))
    return [ r[0] for r in all_in_list], np.array([ r[1:] for r in all_in_list])
    
def sub_2(inum, seuil=0.9):
    fname = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/output/domain_cais_1_2_len%d.lst"%inum
    (label_list, orf1_list, orf2_list) = load_domain_orf1_orf2_list(fname)
    domain_orf1_orf2_list = filter_gcais_less_important(label_list, orf1_list, orf2_list, seuil)

inum=int(sys.argv[1])
    
sub_1(inum)
# sub_2(inum)
