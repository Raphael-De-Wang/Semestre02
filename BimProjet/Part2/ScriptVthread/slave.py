#!env python

import os
import rpyc 
import socket
import pickle
import operator
import subprocess
import numpy as np
import os.path as op
from Bio.Blast.Applications import NcbiblastnCommandline

from BioDist import portNum

def sortDictByValue(dic):
        return sorted(dic.items(), key=operator.itemgetter(1),reverse=True)

def rstFormat(line):
    # Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    qid, sid, idenPercet, alignLen, mismatches, gap, qStart, qEnd, sStart, sEnd, evalue, bitScore = line.split()    
    return ( qid, sid, idenPercet, alignLen, mismatches, gap, qStart, qEnd, sStart, sEnd, evalue, bitScore )

class BioService(rpyc.Service):
    def on_connect(self):
        self.nivExprDict = {}
        print "connected..."
        
    def on_disconnect(self):
        print "disconnected..."
        
    def exposed_slaveHostname(self):
        return socket.gethostname()
    
    def exposed_isFileExist(self,fpath):
        if op.exists(fpath):
            return op.isfile(fname)
        else:
            return False
        
    def exposed_blastHits(self,dbName,qryName,rstName,volPath):
        '''
        blastn [-h] [-help] [-import_search_strategy filename]
        [-export_search_strategy filename] [-task task_name] [-db database_name]
        [-dbsize num_letters] [-gilist filename] [-seqidlist filename]
        [-negative_gilist filename] [-entrez_query entrez_query]
        [-db_soft_mask filtering_algorithm] [-db_hard_mask filtering_algorithm]
        [-subject subject_input_file] [-subject_loc range] [-query input_file]
        [-out output_file] [-evalue evalue] [-word_size int_value]
        [-gapopen open_penalty] [-gapextend extend_penalty]
        [-perc_identity float_value] [-qcov_hsp_perc float_value]
        [-xdrop_ungap float_value] [-xdrop_gap float_value]
        [-xdrop_gap_final float_value] [-searchsp int_value] [-max_hsps int_value]
        [-sum_stats bool_value] [-penalty penalty] [-reward reward] [-no_greedy]
        [-min_raw_gapped_score int_value] [-template_type type]
        [-template_length int_value] [-dust DUST_options]
        [-filtering_db filtering_database]
        [-window_masker_taxid window_masker_taxid]
        [-window_masker_db window_masker_db] [-soft_masking soft_masking]
        [-ungapped] [-culling_limit int_value] [-best_hit_overhang float_value]
        [-best_hit_score_edge float_value] [-window_size int_value]
        [-off_diagonal_range int_value] [-use_index boolean] [-index_name string]
        [-lcase_masking] [-query_loc range] [-strand strand] [-parse_deflines]
        [-outfmt format] [-show_gis] [-num_descriptions int_value]
        [-num_alignments int_value] [-line_length line_length] [-html]
        [-max_target_seqs num_sequences] [-num_threads int_value] [-remote]
        [-version]
        '''
        # blastn -db RNA_Marker113 -query DNA.fa -outfmt 6 -out output_file
        # subprocess.check_output(["blastn", "-db", dbName, "-query", qryName, "-outfmt", "6", "-out", rstName])
        cline = NcbiblastnCommandline(query=qryName, db=dbName, out=rstName, outfmt=6)
        os.chdir(volPath)
        return os.system(str(cline))
    
    def exposed_hitsToNivExpr(self,rstName):
        handle = open(rstName, 'r')
        for line in handle:
            qid = rstFormat(line)[0]
            if self.nivExprDict.has_key(qid):
                self.nivExprDict[qid] += 1
            else:
                self.nivExprDict[qid]  = 1
        handle.close()
        
    def exposed_sortNivExprDict(self):
        self.indice = 0
        self.sortedNivExprList = sortDictByValue(self.nivExprDict)

    def exposed_indiceAddOne(self):
        self.indice += 1

    def exposed_getNivExpr(self):
        if self.indice < len(self.sortedNivExprList):
            return self.sortedNivExprList[self.indice]
        else:
            return ('',0)
                
if __name__ == "__main__":
    from rpyc.utils.server import ThreadedServer
    t = ThreadedServer(BioService, port = portNum)
    print "Starting..."
    t.start()
