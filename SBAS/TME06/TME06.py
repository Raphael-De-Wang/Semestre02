#!env python
from TME06Lib import *

# read Fasta
# fname = "astral-scopedom-seqres-sel-gs-bib-90-2.05.fa"
fname = "astral-scopedom-seqres-sel-gs-bib-30-2.05.fa"
fb = open(fname)
seqDict = readFasta(fb)
fb.close()

# group sequences in family
famDict = groupSeqFamily(seqDict)

# remove family less than certain population
# seuil = 30
# filterSeqFamily(famDict,seuil)
# print "[%d] Families have more than [%d] members\n"%(len(famDict),seuil)

seuil = 50
filterSeqFamily(famDict,seuil)
print "[%d] Families have more than [%d] members\n"%(len(famDict),seuil)

print "Family Names: ", famDict.keys()

# write fasta in family
famFastaPath = "famFasta/"
alignmentPath= "famAlignment/"
hmmPath      = "famHMM/"
outputPath   = "famSearchOutput/"
searchPath   = "famSearch/"

famTrainTestSize = {}

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
        famTrainTestSize[key] = [len(seqListTrain),len(seqListTest)]
        # clustalW alignment
        clustalW(key,famFastaPath,alignmentPath)
        # hmm build -> .hmm model
        hmmBuild(key,alignmentPath,hmmPath)
    save_dict(famTrainTestSize,'famTrainTestSize')
    
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
    famTrainTestSize = load_dict('famTrainTestSize')
    TP,FP = 0,0,0
    testsize = 0
    for key in famDict:
        print key
        testsize += famTrainTestSize.get(key)[1]
        for testFamily in famDict.keys():
            print "\t",testFamily
            paser = parseHmmer3Tab(key,testFamily,searchPath)
            for query in paser:
                # print "\t\t",query.hits[0].bitscore
                print "\t\t",len(query.hits),"hits"
                hits = len(query.hits)
                if key == testFamily:
                    TP += hits
                else:
                    FP += hits
    if TP + FP <> 0:
        FN = testsize - TP - FP
        prec = precision(TP,FP)
        reca = recall(TP,FN)
        fsco = fScore(TP,FP,FN)
        print "family[%s]: precision [%f], recall [%f], fScore [%f]"%(key,prec,reca,fsco)

def RocCurv():
    # step 3 : analyse
    # -- true positive (TP), eqv. with hit
    # -- true negative (TN), eqv. with correct rejection
    # -- false positive (FP), eqv. with false alarm, Type I error
    # -- false negative (FN), eqv. with miss, Type II error
    famTrainTestSize = load_dict('famTrainTestSize')
    testsize = 0
    TP,FP = [],[]
    for key in famDict:
        for testFamily in famDict.keys():
            paser = parseHmmer3Tab(key,testFamily,searchPath)
            for query in paser:
                for hit in query:
                    if hit.query_id == testFamily:
                        TP.append(hit.bitscore)
                    else:
                        FP.append(hit.bitscore)
        TP.sort(reverse=True)
        FP.sort(reverse=True)
        plotRocCurv(key,np.array(TP),np.array(FP),famTrainTestSize.get(key)[1],10,'RocCurv')

step1()
step2()
step3()
RocCurv()
