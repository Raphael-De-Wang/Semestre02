#!env python

import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

seuil = int(sys.argv[1])

fname = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/AT_arc_metatrans.filtered.fasta.cleanup'
oname = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/AT_arc_metatrans.filtered.fasta.cleanup.len'+str(seuil)

ihandle = open(fname)
ohandle = open(oname,'w')

for record in SeqIO.parse(ihandle, 'fasta'):
    if len(record.seq) > seuil * 3 :
        SeqIO.write([record], ohandle, "fasta")
        
ohandle.close()
ihandle.close()
