#!env python

import pylab
import pickle
import operator
import subprocess
import numpy as np
import os.path as op
import matplotlib as mpl
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Align.Applications import ClustalwCommandline

def readFasta(fb):
    seqDict = {}
    for record in SeqIO.parse(fb, "fasta") :
        desc = record.description
        record.desc = descriptionDict(desc)
        if not seqDict.has_key(record.id):
            seqDict[record.id] = record
        else:
            raise ValueError("Duplicated Sequence ID [%d] in fasta. "%(record.id))
    return seqDict


