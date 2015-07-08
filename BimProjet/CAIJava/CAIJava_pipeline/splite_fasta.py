#!env python

from Bio import SeqIO

tgt = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Marine/Chlorophyta/Chlorella_vulgaris.fsa_nt' #  Picochlorum.fsa_nt  Trebouxia_gelatinosa.fsa_nt
# fname_base = '/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Marine/Chlorophyta/Chlorella_vulgaris_%02d.fsa_nt'
fname_base = '/tmp/Chlorella_vulgaris_%02d.fsa_nt'

# tgt = '../source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum_CDS.fasta'
# fname_base = '../source/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum_CDS_sep_%02d.fasta'

handle = open(tgt)
data = SeqIO.parse(handle,"fasta")

i = 0
sep_count = 1
record_list = []
for record in data:
    if i > 100:
        output_handle = open(fname_base%sep_count,'w')
        SeqIO.write(record_list,output_handle,'fasta')
        output_handle.close()
        output_handle = None
        record_list = []
        sep_count += 1
        i = 0
    record_list.append(record)
    i += 1

if i > 0 :
    output_handle = open(fname_base%sep_count,'w')
    SeqIO.write(record_list,output_handle,'fasta')
    output_handle.close()
    output_handle = None

