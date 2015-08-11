#!env python

import sys
import csv
import math
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def CalculateCAIValue(seq, wvalues):
    gcai_log = 0
    L    = 1.0 * len(seq)/3
    miss = 0
    for ind in range(0,len(seq),3):
        key = str(seq.lower()[ind:ind+3])
        try:
            val = wvalues[key]
        except:
            miss += 1.
            # print "Warning: Encounter Key Inconu. [%s] "%key
            continue
        if val > 0:
            gcai_log += math.log(val)
        else:
            gcai_log += math.log(0.01)
    # if miss > 20:
    #     raise ValueError("Too Many Missed Keys.")
    if L == miss :
        return -1
    return math.exp(gcai_log/(L-miss))

def wvalue_to_dict(wvalue_handler):
    wv_dict = {}
    reader = csv.reader(wvalue_handler)
    for [key,value] in reader :
        if not wv_dict.has_key(key):
            wv_dict[key] = float(value)
        else:
            raise ValueError("Weight Key and Value Conflict.")
    return wv_dict

domFasta_handler = open(sys.argv[1]) # fasta file path
wvalue_handler   = open(sys.argv[2]) # wvalue file path
wvalue_dict = wvalue_to_dict(wvalue_handler)
wvalue_handler.close()

gcai_handler = open(sys.argv[3], 'w') # gcai output file path
csv_writer   = csv.writer(gcai_handler, delimiter='\t')

gcai_list = []

for record in SeqIO.parse(domFasta_handler, "fasta") :
    gcai = CalculateCAIValue(record.seq, wvalue_dict)
    if gcai >= 0 :
        gcai_list.append([record.name,gcai])

gcai_list = sorted(gcai_list, key=operator.itemgetter(1),reverse=True)
    
for record in gcai_list:    
    csv_writer.writerow(record)

domFasta_handler.close()
gcai_handler.close()    

