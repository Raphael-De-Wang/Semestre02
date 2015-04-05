#!env python
import pylab
import operator
import subprocess
import numpy as np
import os.path as op
import matplotlib as mpl
import matplotlib.pyplot as plt

from Bio import SeqIO
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
    
def aligClustalW(familyName,inPath,outPath):
    clustalwCMD = subprocess.check_output(["which","clustalw"])
    clustalw_cline = ClustalwCommandline(clustalwCMD, infile=trainFastaFileName(inPath,familyName), outfile=op.join(outPath,familyName+'.aln'))
    clustalw_cline()

