#!env python
import sys
from Bio import SeqIO,Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

fin  = sys.argv[1]
fout = sys.argv[2]

in_handle  = open(fin)
out_handle = open(fout,'w')
data = SeqIO.parse(in_handle,'fasta')
# trans_list = []

for record in data:
    trans = SeqRecord(Seq.translate(record.seq[0:]),id=record.id,name=record.name,description=record.description)
    SeqIO.write([trans], out_handle, "fasta")
    # trans_list.append(trans)
    
out_handle.close()

