#!env python

import sys
from Bio import SeqIO

# tgt = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum.gb'
# tgt = '/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/source/Marine/Bacillariophyta/Thalassiosira_pseudonana.gb'

input  = sys.argv[1]
output = sys.argv[2]

'''
for idx,record in enumerate(SeqIO.parse(open(tgt), "genbank")) :
    i = 0
    for f in record.features : 
        if f.type == 'CDS' :
            i += 1
    print  idx,record.annotations['accessions'][0],' has %d CDS '%i
'''

i = 0
seq_list = []
input_handle = open(input)
for idx,record in enumerate(SeqIO.parse(input_handle, "genbank")) :
    for f in record.features :
        if f.type <> 'CDS' :
            continue
        seq = f.extract(record)
        if len(record) > len(seq) :
            seq.dbxrefs     = f.qualifiers['db_xref']
            seq.description = ''# f.qualifiers['product'][0]
            # seq.name = f.qualifiers['gene']
            try:
                seq.name = '|'.join(f.qualifiers['db_xref'] + f.qualifiers['locus_tag'] + f.qualifiers['protein_id'] + f.qualifiers['product']).replace(':','|')
            except:
                try:
                    seq.name = '|'.join(f.qualifiers['db_xref'] + f.qualifiers['locus_tag']).replace(':','|')
                except:
                    seq.name = '|'.join(f.qualifiers['db_xref']).replace(':','|')
            seq.id   = seq.name
            seq_list.append(seq)
            
# output_handle = open('/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum.fasta','w')
# output_handle = open(tgt.replace('.gb','_CDS.fasta'),'w')

output_handle = open(output,'w')
SeqIO.write(seq_list, output_handle, "fasta")
    
input_handle.close()
output_handle.close()
