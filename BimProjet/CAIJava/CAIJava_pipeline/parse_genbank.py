#!env python

import numpy as np
from Bio import SeqIO
from subprocess import call

workpath = "/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/"

'''
for seq_record in SeqIO.parse(workpath+"CAIJava/source/Bacillariophyta_full.gb", "genbank"):
    handle = open('/tmp/'+seq_record.id+'.gbk','w')
    SeqIO.write(seq_record,handle,'genbank')
    handle.close()
'''

# src_handle = open(workpath+"CAIJava/source/Bacillariophyta_full.gb")
# src_handle = open(workpath+"CAIJava/source/Bacillariophyta.gb")
src_handle = open('/home/raphael/Documents/UPMC_BIM/Semestre02/TME/BimProjet/CAIJava/source/Thalassiosira_pseudonana/Thalassiosira_pseudonana.gb')

flag = True
out_handle = None
source = None
accession  = None

for line in src_handle:
    if 'SOURCE' in line:
        source = '_'.join(line.strip().split()[1:])
    if flag == True and len(line.strip()) <> 0 :
        flag = False
        accession = line.split()[1]
        out_handle = open('tmpSource/'+accession+'.gbk','w')
    if out_handle <> None:
        out_handle.write(line)
    if "//" == line.strip() :
        flag = True
        call(['mv','tmpSource/'+accession+'.gbk', 'tmpSource/'+accession+'_'+source+'.gbk'])
        # print line
        out_handle.close()
        out_handle = None
    '''
    if source <> None and accession <> None:
        print source,',',accession
        source = None
        accession  = None
    '''
