#!env python

from Bio import SeqIO,Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
# 'source/Marine/Chlorophyta/Chlorella_vulgaris.fsa_nt',
# 'source/Marine/Chlorophyta/Picochlorum.fsa_nt',
# 'source/Marine/Chlorophyta/Trebouxia_gelatinosa.fsa_nt',

tgt_list = ['source/Marine/Dinoflagellata/Symbiodinium_minutum.fsa_nt',
            'source/Marine/Rhodophyta/Porphyridium_purpureum.fsa_nt',
            'source/Marine/Streptophyta/Klebsormidium_flaccidum.fsa_nt',
            'source/Marine/Streptophyta/Paramecium_biaurelia.fsa_nt',
            'source/Marine/Streptophyta/Paramecium_caudatum.fsa_nt',
            'source/Marine/Streptophyta/Paramecium_sexaurelia.fsa_nt',
            'source/Marine/Streptophyta/Picochlorum.fsa_nt',
            'source/Marine/Streptophyta/Spirodela_polyrhiza.fsa_nt',
            'source/Marine/Streptophyta/Tetrahymena_borealis.fas_nt',
            'source/Marine/Streptophyta/Tetrahymena_elliotti.fsa_nt',
            'source/Marine/Streptophyta/Tetrahymena_malaccensis.fsa_nt']

path_base = '/users/Etu9/3404759/Workspace/Semestre02/BimProjet/CAIJava/'

for tgt in tgt_list:
    input_handle  = open(path_base+tgt)
    output_fname  = path_base+tgt.split('.')[0]+'_translate.fasta'
    print output_fname
    output_handle = open(output_fname,'w')
    data = SeqIO.parse(input_handle,'fasta')
    for record in data:
        frame01 = SeqRecord(Seq.translate(record.seq[0:]),id=record.id+'rdf1|',
                            name=record.name+'rdf1|',description=record.description)
        frame02 = SeqRecord(Seq.translate(record.seq[1:]),id=record.id+'rdf2|',
                            name=record.name+'rdf2|',description=record.description)
        frame03 = SeqRecord(Seq.translate(record.seq[2:]),id=record.id+'rdf3|',
                            name=record.name+'rdf3|',description=record.description)
        frame04 = SeqRecord(Seq.translate(record.reverse_complement().seq[0:]),id=record.id+'rdf4|',
                            name=record.name+'rdf4|',description=record.description)
        frame05 = SeqRecord(Seq.translate(record.reverse_complement().seq[1:]),id=record.id+'rdf5|',
                            name=record.name+'rdf5|',description=record.description)
        frame06 = SeqRecord(Seq.translate(record.reverse_complement().seq[2:]),id=record.id+'rdf6|',
                            name=record.name+'rdf6|',description=record.description)
        SeqIO.write([frame01,frame02,frame03,frame04,frame05,frame06], output_handle, "fasta")
    output_handle.close()

