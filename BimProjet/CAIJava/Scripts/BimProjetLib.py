#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt
import pylab

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import FastaWriter

WORKPATH="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/"

#### Format ####
def toStringGNameRFBeginEndDomain(gName,RF,begin,end,domain):
    return "%s_%s_%d_%d_%s"%(gName,RF,begin,end,domain)

def gNameRFBeginEndDomainToStrArray(ID):
    (gName,RF,begin,end,domain) = ID.split("_")
    return (gName,RF,int(begin),int(end),domain)
    
def Step2FormatSepBestDomain(line):
    score = float(line[0])
    gId   = line[4]
    gName,ARC,RF = gId.split("__")
    domain = line[5]
    RF = RF.split("_")[-1]
    reverse,begin,end = indiceBeginEnd(RF,int(line[2]),int(line[3]))
    gId   = toStringGNameRFBeginEndDomain(gName,RF,begin,end,domain)
    return (score,gName,domain,gId,ARC,RF,reverse,begin,end,"")

def Step2FormatSepArchs(line):
    score = float(line[0])
    gId   = line[3]
    domain= line[4]
    gName,ARC,RF = gId.split("__")
    RF = RF.split("_")[-1]
    reverse,begin,end = indiceBeginEnd(RF,int(line[1]),int(line[2]))
    gId   = toStringGNameRFBeginEndDomain(gName,RF,begin,end,domain)
    return (score,gName,domain,gId,ARC,RF,reverse,begin,end,"")

def readFASTA(fname="AT_arc_metatrans.filtered.fasta"):
    return SeqIO.index(fname, "fasta")

def indiceBeginEnd(readFrame,begin,end):
    rf = int(readFrame[-1])
    return (rf>3,(begin-1)*3+rf-1,end*3+rf-1)

def cleanUpFasta(fname,fastaDict,step2List,Step2FormatSepFunc,seuil=1e-3):
    with open(fname,'w') as fb:
        writer = FastaWriter(fb)
        writer.write_header()
        for line in step2List:
            try:
                score,gName,domain,gId,ARC,RF,reverse,begin,end,desc = Step2FormatSepFunc(line)
                if score > seuil:
                    # print "[%s] score [%f] > seuil [%f].\n"%(gName,score,seuil)
                    continue
                code = fastaDict[gName].seq.tostring()
                if reverse:
                    code = code[::-1]
                record = SeqRecord(Seq(code[begin:end],generic_dna),name=gName,id=gId,description=desc)
                writer.write_record(record)
            except KeyError:
                print "[%s] not exists in fasta dictionary.\n"%gName
                continue
        writer.write_footer()
    
def readgCAIsKeyGName(fname=WORKPATH+"output/cais.lst"):
    fb = open(fname)
    CAI_dict = dict()
    for line in fb:
        line_array = line.strip().split("\t")
        gName      = line_array[0].split('_')[0]
        gCAI       = line_array[1]
        try:
            if CAI_dict.get(gName) <> None:
                value = CAI_dict[gName]
                CAI_dict[gName].append(float(gCAI))
            else:
                CAI_dict[gName] = [float(gCAI)]
        except:
            if gCAI == "CAI":
                continue
            raise
    fb.close()
    return CAI_dict

def readgCAIs(fname=WORKPATH+"output/cais.lst"):
    fb = open(fname)
    CAI_dict = dict()
    for line in fb:
        line_array = line.strip().split("\t")
        gName      = line_array[0]
        gCAI       = line_array[1]
        try:
            CAI_dict[gName] = float(gCAI)
        except:
            if gCAI == "CAI":
                continue
            raise
    fb.close()
    return CAI_dict

def gCAIsDictToDomDict(CAI_dict,domIdDict,keys):
    gCAIsDomDict = dict()
    for key in keys:
        ids = domIdDict[key]
        gCAIsDomDict[key] = []
        for i in ids :
            gCAIsDomDict[key].append(CAI_dict.get(i))
    return gCAIsDomDict
    
def readClstr(fname):
    fb       = open(fname)
    CLS_dict = dict()
    nivExp   = 0
    gname    = None
    for line in fb:
        if line[0] == '>':
            if gname == None:
                continue
            CLS_dict[gname] = nivExp
            nivExp = 0
            gname  = None
            continue
        larray = line.strip().split()
        nivExp = int(larray[0])+1
        if larray[-1] == '*':
            gname = larray[2].replace('>','').replace('...','')
    if gname<>None and nivExp<>0:
        CLS_dict[gname] = nivExp
    fb.close()
    return CLS_dict

def read2step(fname):
    fb = open(fname)
    M = np.array([ line.strip().split() for line in fb ])
    fb.close()
    return M

def domainIdDict(idList):
    did = dict()
    for ID in idList:
        gName,RF,begin,end,domain = gNameRFBeginEndDomainToStrArray(ID)
        if did.has_key(domain):
            if ID not in did.get(domain):
                did[domain].append(ID)
        else:
            did[domain] = [ID]
    return did

def domainGenomeDict(idList):
    dgs = dict()
    for ID in idList:
        gName,RF,begin,end,domain = gNameRFBeginEndDomainToStrArray(ID)
        if dgs.has_key(domain):
            if gName not in dgs.get(domain):
                dgs[domain].append(gName)
        else:
            dgs[domain] = [gName]
    return dgs

def accumulateNivExprDom(nivExprDomDict,domIdDict):
    acc = dict()
    for key in domIdDict:
        acc[key] = 0
        for gName in domIdDict[key]:
            acc[key] = acc[key] + nivExprDomDict[gName]
    return acc

def sortDictByValue(dic):
    return sorted(dic.items(), key=operator.itemgetter(1),reverse=True)

def sortListByCol(tList,colNum):
    return sorted(tList,key=lambda lst: lst[colNum])

def binSearch(tList,comFunc,match):
    haut,bas = len(tList)-1,0
    while haut >= bas :
        mid  = int(round((haut+bas)*1./2))
        bias = comFunc(tList[mid],match)
        if bias == 0:
            return mid
        elif bias > 0:
            bas = mid + 1
        else:
            haut = mid - 1

def step2GnameComp(step2Elmt,match):
    try:
        if match in step2Elmt[3]:
            return 0
        elif step2Elmt[3][:len(match)] > match:
            return 1
        else:
            return -1
    except:
        print step2Elmt
        raise
        exit()

def autolabel(ax,rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height), ha='center', va='bottom')

def plotBarChart(plotList,labels,figName=None):
    fig, ax = plt.subplots()
    fig.set_figheight(fig.get_figheight()+2)
    rects   = ax.bar(np.arange(len(plotList)), plotList, color='blue')
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Niveau d\'Expression')
    ax.set_title('Sorted Niveau d\'Expression par domaine')
    # ax.set_xticks(ind+width)
    ax.set_xticklabels(labels,rotation=90)
    # autolabel(ax,rects)
    if figName:
        plt.savefig(figName)
    else:
        plt.show()
    plt.close()
    
def plotgCAIsDomain(X,Y,figName=None):
    pylab.figure(figsize=(20,4))
    pylab.hist2d(Y,X,bins=[len(np.unique(Y)),50])
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.colorbar()
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close()
    
def plotgCAIsDomain2(X,Y,figName=None):
    T = np.arctan(X)#,np.ones(len(X)))
    pylab.figure(figsize=(20,4))
    pylab.scatter(Y,X, s=75, c=T, alpha=.5)
    pylab.ylim(0,1.2), pylab.yticks([])
    pylab.xlim(0,len(np.unique(Y))), pylab.xticks([])
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close()
    
def plotgCAISNivExpr(gCAIs,nivExpr,figName=None):
    fig = plt.figure()
    # plt.subplot(223)
    plt.ylim(ymax=1.2)
    # plt.xlabel("Genome")
    # plt.ylabel("Score")
    plt.plot(gCAIs,'ro',label="gCAIs")
    plt.plot(nivExpr*1./max(nivExpr),label="Niveau d'Expression")
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.legend(loc=2, borderaxespad=0.)
    if figName:
        plt.savefig(figName)
    else:
        plt.show()
    plt.close()
    
def plotHist(scores,figName=None):
    fig = plt.figure()
    plt.hist(scores)
    if figName:
        plt.savefig(figName)
    else:
        plt.show()
    plt.close()
        
def saveIDgCAInivExpr(fname,data):
    fb = open(fname,'w')
    # fb.write("ID\tgCAIs\tNivExpr\n")
    for line in data:
        fb.write("%s\t%s\t%s\n"%(line[0],str(line[1]),str(line[2])))

def saveIDgCAInivExprMarked(fname,data):
    fb = open(fname,'w')
    # fb.write("ID\tgCAIs\tNivExpr\tMarked\n")
    for line in data:
        fb.write("%s\t%s\t%s\t%s\n"%(line[0],str(line[1]),str(line[2]),str(line[3])))

def mergegCAIsNivExprByGname(sortList, mergeDict):
    return np.array([ np.append(line,mergeDict.get(line[0])) for line in sortList])

def mergeStep2ByGname(sortList,step2M):
    return np.array([ np.append(line,binSearch(step2M,step2GnameComp,line[0])>=0) for line in sortList])

def nivExprDomain(step2MList, colNum):
    nivExprDict = dict()
    for record in step2MList:
        domain = record[colNum]
        ne = nivExprDict.get(domain) 
        if ne is None:
            nivExprDict[domain] = 1
        else:
            nivExprDict[domain] = ne + 1
    return nivExprDict

def getGIDList(step2MList):
    GIDs = []
    for line in step2MList:
        (score,gName,domain,gId,ARC,RF,reverse,begin,end,desc) = Step2FormatSepArchs(line)
        GIDs.append(toStringGNameRFBeginEndDomain(gName,RF,begin,end,domain))
    return np.array(GIDs)

#### trier la donnee ####
def trierFastaByDomain(tgtDomain,fastaDict,step2List,writeFileName,formatFunc):
    fb     = open(writeFileName,'w')
    writer = FastaWriter(fb)
    writer.write_header()
    for record in step2List:
        score,gName,domain,gID,ARC,RF,reverse,begin,end,desc = formatFunc(record)
        if domain == tgtDomain:
            if fastaDict.get(gID) <> None:
                writer.write_record(fastaDict.get(gID))
            '''
            else:
                print "[%s] n'existe pas dans le fiche"%(gID)
            '''
    writer.write_footer()
    fb.close()

