#!env python
from TME06Lib import *

# read Fasta
fname = "astral-scopedom-seqres-sel-gs-bib-90-2.05.fa"
fb = open(fname)
seqDict = readFasta(fb)
fb.close()

# group sequences in family
famDict = groupSeqFamily(seqDict)

# remove family less than 30 population
seuil = 100
filterSeqFamily(famDict,seuil)
print "[%d] Families have more than [%d] members\n"%(len(famDict),seuil)

# write fasta in family
famFastaPath = "famFasta/"
alignmentPath= "famAlignment/"

for key in famDict:
    seqList = famDict.get(key)
    seqListTrain,seqListTest = sepEchantillon(seqList,0.75)
    fbTrain = open(trainFastaFileName(famFastaPath,key),'w')
    fbTest  = open(testFastaFileName(famFastaPath,key),'w')
    writeFasta(fbTrain,seqListTrain)
    writeFasta(fbTest ,seqListTest)
    fbTrain.close()
    fbTest.close()
    # clustalW
    aligClustalW(key,famFastaPath,alignmentPath)
    break

