#!env python

from Bio import SeqIO

# tgt = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum.gb'
tgt = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Thalassiosira_pseudonana/Thalassiosira_pseudonana.gb'

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
input_handle = open(tgt)
for idx,record in enumerate(SeqIO.parse(input_handle, "genbank")) :
    for f in record.features :
        if f.type <> 'CDS' :
            continue
        seq = f.extract(record)
        if len(record) > len(seq) :
            seq.dbxrefs     = f.qualifiers['db_xref']
            seq.description = ''# f.qualifiers['product'][0]
            # seq.name = f.qualifiers['gene']
            seq.name = f.qualifiers['locus_tag'][0]
            seq.id   = f.qualifiers['db_xref'][0]
            seq_list.append(seq)
            
# output_handle = open('/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum.fasta','w')
output_handle = open('/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Thalassiosira_pseudonana/Thalassiosira_pseudonana.fasta','w')
SeqIO.write(seq_list, output_handle, "fasta")
    
input_handle.close()
output_handle.close()
