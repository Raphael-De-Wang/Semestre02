#!env python
from TME06Lib import *

# read Fasta
fname = "astral-scopedom-seqres-sel-gs-bib-90-2.05.fa"
fb = open(fname)
seqDict = readFasta(fb)
fb.close()

# group sequences in family
famDict = groupSeqFamily(seqDict)

# remove family less than 100 population
seuil = 150
filterSeqFamily(famDict,seuil)
print "[%d] Families have more than [%d] members\n"%(len(famDict),seuil)

# write fasta in family
famFastaPath = "famFasta/"
alignmentPath= "famAlignment/"
hmmPath      = "famHMM/"
outputPath   = "famSearchOutput/"
searchPath   = "famSearch/"

def step1():
    for key in famDict:
        seqList = famDict.get(key)
        # step 1 : seperate training set and testing set
        seqListTrain,seqListTest = sepEchantillon(seqList,0.75)
        if len(seqListTrain) < 1 or len(seqListTest) < 1:
            continue
        fbTrain = open(trainFastaFileName(famFastaPath,key),'w')
        writeFasta(fbTrain,seqListTrain)
        fbTrain.close()
        fbTest  = open(testFastaFileName(famFastaPath,key),'w')
        writeFasta(fbTest ,seqListTest)
        fbTest.close()
        # clustalW alignment
        aligClustalW(key,famFastaPath,alignmentPath)
        # hmm build -> .hmm model
        hmmBuild(key,alignmentPath,hmmPath)
    
def step2():
    # step 2 : hmmsearch
    for key in famDict:
        for testFamily in famDict.keys():
            hmmSearch(key,hmmPath,testFamily,famFastaPath,outputPath,searchPath)
    
def step3():
    # step 3 : analyse
    # -- true positive (TP), eqv. with hit
    # -- true negative (TN), eqv. with correct rejection
    # -- false positive (FP), eqv. with false alarm, Type I error
    # -- false negative (FN), eqv. with miss, Type II error
    for key in famDict:
        for testFamily in famDict.keys():
            paser = parseHmmer3Tab(key,testFamily,searchPath)
            for query in paser:
                print query.hits[0].bitscore

# step1()
# step2()
step3()
    
