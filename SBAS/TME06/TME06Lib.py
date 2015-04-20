#!env python
import pylab
import pickle
import operator
import subprocess
import numpy as np
import os.path as op
import matplotlib as mpl
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Align.Applications import ClustalwCommandline


def descriptionDict(desc):
    dDict = {}
    dList = desc.split()
    dDict["id"]      = dList[0]
    dDict["family"] = dList[1]
    return dDict
    
def readFasta(fb):
    seqDict = {}
    for record in SeqIO.parse(fb, "fasta") :
        desc = record.description
        record.desc = descriptionDict(desc)
        if not seqDict.has_key(record.id):
            seqDict[record.id] = record
        else:
            raise ValueError("Duplicated Sequence ID [%d] in fasta. "%(record.id))
    return seqDict

def writeFasta(fb,seqList):
    if len(seqList) <= 0:
        raise ValueError("No data to Persist.")
    writer = FastaWriter(fb)
    writer.write_header()
    for record in seqList:
        writer.write_record(record)
    writer.write_footer()
    
def groupSeqFamily(seqDict):
    famDict = {}
    for key in seqDict:
        record = seqDict.get(key)
        fam    = record.desc.get("family")
        if famDict.has_key(fam):
            famDict.get(fam).append(record)
        else:
            famDict[fam] = []
    return famDict

def filterSeqFamily(famDict,seuil=30):
    rList = []
    for key in famDict:
        if len(famDict.get(key)) < seuil:
            rList.append(key)
    for key in rList:
        famDict.pop(key)

def sepEchantillon(seqList,trainRate):
    indice = range(len(seqList))
    trainNum = int(trainRate*len(seqList))
    np.random.shuffle(seqList)
    return ( seqList[:trainNum], seqList[trainNum:] )

def trainFastaFileName(path,familyName):
    return op.join(path,familyName+'.train.fa')
    
def testFastaFileName(path,familyName):
    return op.join(path,familyName+'.test.fa')

def alnFastaFileName(path,familyName):
    return op.join(path,familyName+'.aln')
    
def hmmModelName(path,familyName):
    return op.join(path,familyName+'.hmm')
    
def outputSearchFileName(path,modelFamilyName,dbFamilyName):
    return op.join(path,modelFamilyName+'-'+dbFamilyName+'.output')
    
def hmmer3TabName(path,modelFamilyName,dbFamilyName):
    return op.join(path,modelFamilyName+'-'+dbFamilyName+'.hmmer3-tab')
    
def clustalW(familyName,inPath,outPath):
    clustalwCMD = subprocess.check_output(["which","clustalw2"]).strip()
    clustalw_cline = ClustalwCommandline(clustalwCMD, infile=trainFastaFileName(inPath,familyName), outfile=alnFastaFileName(outPath,familyName))
    return subprocess.check_output(str(clustalw_cline).split())

def hmmBuild(familyName,inPath,outPath):
    hmmBuildCMD = subprocess.check_output(["which","hmmbuild"]).strip()
    rst = subprocess.check_output([hmmBuildCMD,hmmModelName(outPath,familyName),alnFastaFileName(inPath,familyName)])
    
def hmmSearch(modelFamilyName,hmmPath,dbFamilyName,famFastaPath,outputPath,searchPath):
    hmmSearchCMD = subprocess.check_output(["which","hmmsearch"]).strip()
    rst = subprocess.check_output([hmmSearchCMD,"-o", outputSearchFileName(outputPath,modelFamilyName,dbFamilyName),
                                    '--tblout', hmmer3TabName(searchPath,modelFamilyName,dbFamilyName),
                                    hmmModelName(hmmPath,modelFamilyName),
                                    testFastaFileName(famFastaPath,dbFamilyName)])

def parseHmmer3Tab(modelFamilyName,dbFamilyName,searchPath):
    return SearchIO.parse(hmmer3TabName(searchPath,modelFamilyName,dbFamilyName),'hmmer3-tab')


def precision(TP,FP):
    return 1.*TP/(TP+FP)

def recall(TP,FN):
    return 1.*TP/(TP+FN)

def fScore(TP,FP,FN):
    return 2.*precision(TP,FP)*recall(TP,FN)/(precision(TP,FP)+recall(TP,FN))

def save_dict(d,fname):
    fb = open(fname,'w')
    pickle.dump(d,fb)
    fb.close()

def load_dict(fname):
    fb = open(fname)
    b = pickle.load(fb)
    fb.close()
    return b

def plotRocCurv(key,TP,FP,P,interv,fname=None):
    vmax = max(np.ravel(TP.tolist(),FP.tolist()))
    x,y = [],[]
    # b110_tp,b110_fp,b110_fn,b110_sc = [],[],[],[]
    ran  = int(vmax/interv+1)
    for i in range(1,ran):
        vtp = sum(TP>i*interv)
        vfp = sum(FP>i*interv)
        vfn = P - vtp - vfp
        if vtp + vfp <> 0:
            prec = precision(vtp,vfp)
            reca = recall(vtp,vfn)
            x.append(reca)
            y.append(prec)
            '''
            if key == 'b.1.1.0':
                print "family[%s] %d: precision [%f], recall [%f]"%(key,i,prec,reca)
                print "TP: %f, FP: %f, FN: %f, Score: %d"%(vtp,vfp,vfn,i*interv)
                b110_tp.append(vtp)
                b110_fp.append(vfp)
                b110_fn.append(vfn)
                b110_sc.append(i*interv)
            '''
    fig = plt.figure()
    plt.title("Roc Curv: %s family"%key)
    plt.xlabel("Specificity(recall %)")
    plt.ylabel("Sensitivity(precision %)")
    plt.ylim((0,1.1))
    plt.plot(x,y)
    if fname <> None:
        plt.savefig(fname.replace('.','-')+'-'+key+'.png')
    else:
        plt.show()
    plt.close(fig)
    '''
    if key == 'b.1.1.0':
        fig = plt.figure()
        indice = b110_sc
        plt.title("TP,FP,FN: %s family"%key)
        plt.xlabel("Score")
        plt.ylabel("Value")
        plt.plot(indice,b110_tp,label='TP')
        plt.plot(indice,b110_fp,label='FP')
        plt.plot(indice,b110_fn,label='FN')
        plt.legend(loc=1)
        plt.savefig(key.replace('.','-')+'TP-FP-FN'+'.png')
        plt.close(fig)
    '''
