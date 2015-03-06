#!/usr/bin/env python
import operator
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import FastaWriter

WORKPATH="/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/"

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
                score,gName,gId,ARC,RF,reverse,begin,end,desc = Step2FormatSepFunc(line)
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
    
def readgCAIs(fname=WORKPATH+"output/cais.lst"):
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
        
def plotgCAISNivExpr(gCAIs,nivExpr,figName=None):
    fig = plt.figure()
    # plt.subplot(223)
    plt.ylim(ymax=1.2)
    plt.xlabel("Genome")
    plt.ylabel("Score")
    plt.plot(gCAIs,label="gCAIs")
    plt.plot(nivExpr*1./max(nivExpr),label="Niveau d'Expression")
    # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.legend(loc=2, borderaxespad=0.)
    if figName:
        plt.savefig(figName)
    else:
        plt.show()
    
def plotHist(scores,figName=None):
    fig = plt.figure()
    plt.hist(scores)
    if figName:
        plt.savefig(figName)
    else:
        plt.show()
    
def saveIDgCAInivExpr(fname,data):
    fb = open(fname,'w')
    fb.write("ID\tgCAIs\tNivExpr\n")
    for line in data:
        fb.write("%s\t%s\t%s\n"%(line[0],str(line[1]),str(line[2])))

def saveIDgCAInivExprMarked(fname,data):
    fb = open(fname,'w')
    fb.write("ID\tgCAIs\tNivExpr\tMarked\n")
    for line in data:
        fb.write("%s\t%s\t%s\t%s\n"%(line[0],str(line[1]),str(line[2]),str(line[3])))

def mergegCAIsNivExprByGname(sortList, mergeDict):
    return np.array([ np.append(line,mergeDict.get(line[0])) for line in sortList])

def mergeStep2ByGname(sortList,step2M):
    return np.array([ np.append(line,binSearch(step2M,step2GnameComp,line[0])>=0) for line in sortList])

def nivExprDomain(step2List, colNum):
    nivExprDict = dict()
    for record in step2List:
        domain = record[colNum]
        ne = nivExprDict.get(domain) 
        if ne is None:
            nivExprDict[domain] = 1
        else:
            nivExprDict[domain] = ne + 1
    return nivExprDict
