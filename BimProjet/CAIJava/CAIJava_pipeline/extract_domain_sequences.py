#!env python

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

class domItem(object) :
    def __init__(self, record):
        tname, dom_accession, tlen, qname, ful_seq_accession, qlen, evalue, ful_seq_score, ful_seq_bias, item_id, item_num, cEvalue,  iEvalue, dom_score, dom_bias, hmmfrom, hmmto, alifrom, alito, envfrom, envto, acc = record[:22];
        self.tname    = tname
        self.dom_accession = dom_accession
        self.tlen     = int(tlen)
        self.qname    = qname
        self.ful_seq_accession = ful_seq_accession
        self.qlen     = int(qlen)
        self.evalue   = float(evalue)
        self.ful_seq_score = float(ful_seq_score)
        self.ful_seq_bias  = float(ful_seq_bias)
        self.item_id  = int(item_id)
        self.item_num = int(item_num)
        self.cEvalue  = float(cEvalue)
        self.iEvalue  = float(iEvalue)
        self.dom_score= float(dom_score)
        self.dom_bias = float(dom_bias )
        self.hmmfrom  = int(hmmfrom)
        self.hmmto    = int(hmmto)
        self.alifrom  = int(alifrom)
        self.alito    = int(alito)
        self.envfrom  = int(envfrom)
        self.envto    = int(envto)
        self.acc      = float(acc)
        self.desc     = ' '.join(record[22:])

domtblout = open(sys.argv[1])
domFasta  = open(sys.argv[2], 'w')
cds       = None

seqDict = {}
record_list = []

ORF6 = False

if len(sys.argv) > 3:
    if cmp(sys.argv[3],"-orf6") == 0:
        cds  = open(sys.argv[1].replace(".domtblout",".fsa_nt"))
        ORF6 = True
else:
    cds = open(sys.argv[1].replace(".domtblout","_CDS.fasta"))

for record in SeqIO.parse(cds, "fasta") :
    name = '|'.join(record.name.split('|')[:2])
    if seqDict.has_key(name) :
        print record.name
        print seqDict[name].name
    seqDict[name] = record

def normal_extract(domtblout,seqDict,seuil=0):
    for line in domtblout :
        if line[0] == '#':
            continue
        dom = domItem(line.split())
        name = '|'.join(dom.qname.split('|')[:2])
        ref = name + '|' + dom.dom_accession + '|' + "%d.%d"%(dom.item_num, dom.item_id)
        if dom.alito - dom.alifrom <= seuil :
            continue
        seq = Seq(str(seqDict[name].seq)[dom.alifrom*3:dom.alito*3], IUPAC.unambiguous_dna)
        record = SeqRecord(seq, id=ref, name=ref, description=dom.desc)
        record_list.append(record)

def orf6_extract(domtblout,seqDict, seuil=0):
    for line in domtblout :
        if line[0] == '#':
            continue
        dom = domItem(line.split())
        name = '|'.join(dom.qname.split('|')[:2])
        ref = name + '|' + dom.dom_accession + '|' + "%d.%d"%(dom.item_num, dom.item_id)
        if dom.alito - dom.alifrom <= seuil :
            continue
        if   cmp("|ORF1|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].seq)[dom.alifrom*3:dom.alito*3], IUPAC.unambiguous_dna)
        elif cmp("|ORF2|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].seq)[dom.alifrom*3+1:dom.alito*3+1], IUPAC.unambiguous_dna)
        elif cmp("|ORF3|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].seq)[dom.alifrom*3+2:dom.alito*3+2], IUPAC.unambiguous_dna)
        elif cmp("|ORF4|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].reverse_complement().seq)[dom.alifrom*3:dom.alito*3], IUPAC.unambiguous_dna)
        elif cmp("|ORF5|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].reverse_complement().seq)[dom.alifrom*3+1:dom.alito*3+1], IUPAC.unambiguous_dna)
        elif cmp("|ORF6|",dom.qname[-6:]) == 0:
            seq = Seq(str(seqDict[name].reverse_complement().seq)[dom.alifrom*3+2:dom.alito*3+2], IUPAC.unambiguous_dna)
        else:
            raise ValueError("Sequence ORF Inconu.")
        record = SeqRecord(seq, id=ref, name=ref, description=dom.desc)
        record_list.append(record)

seuil = 30
        
if ORF6:
    orf6_extract(domtblout,seqDict,seuil)
else:
    normal_extract(domtblout,seqDict,seuil)

SeqIO.write( record_list, domFasta, "fasta")
domFasta.close()

