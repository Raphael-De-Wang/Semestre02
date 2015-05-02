#!env python

import numpy as np
from Bio import SeqIO

code_dict = {'A':0,'C':1,'G':2,'T':3,
             'a':0,'c':1,'g':2,'t':3}

def readFasta(handle):
    seqList = []
    for record in SeqIO.parse(handle, "fasta") :
        seqList.append(record)
    return seqList

def P(seq,pos,pt,W):
    prob = 1
    for p in range(0,pos):
        prob *= pt[code_dict[seq[p]],0]
    for p in range(pos,pos+W):
        prob *= pt[code_dict[seq[p]],p-pos+1]
    for p in range(pos+W,len(seq)):
        prob *= pt[code_dict[seq[p]],0]
    return prob
        
def likelihood_log(seqList,Z,pt,W):
    ll = 0
    for i,record in enumerate(seqList):
        pos = np.argmax(Z[i])
        # print "Z[i]: \n", Z[i]
        # print "pos : ", pos
        seq = record.seq
        ll += np.log(P(seq,pos,pt,W)/(len(seq)-W+1))
    return ll

def init_p0_random(seqList,W): 
    p0 = np.ones((4,W+1)) / 6
    p0[:,0] = np.ones(4)*0.25
    i = np.random.randint(len(seqList))
    j = np.random.randint(len(seqList[i].seq)-W)
    print "seq%d, position %d"%(i,j)
    print "random motif : ", seqList[i].seq[j:j+W]
    for indice,c in enumerate(seqList[i].seq[j:j+W]):
        p0[code_dict[c],indice+1] = 0.5
    return p0

def E_step(seqList,pt,W):
    Zt = np.array([[ P(record.seq,j,pt,W) for j in range(len(record.seq)-W+1) ] for i,record in enumerate(seqList) ])
    # normalisation de Zt
    for i,s in enumerate(Zt.sum(axis=1)):
        Zt[i] = Zt[i]/s
    return Zt

def Nc(seq,c):
    return sum(np.array(list(seq)) == c)

def M_step(seqList,W,Zt):
    # Calcul de Pck
    pt = np.ones((4,W+1))
    for i,record in enumerate(seqList):
        for pos,c in enumerate(record.seq):
            L = len(record.seq)
            if pos + 1 >= W :
                if pos <= L - W + 1:
                    jRange = range(0,W)
                else:
                    jRange = range(pos+W-L-1,W)
            else :
                jRange = range(0,pos+1)
            for j in jRange:
                pt[code_dict[c],j+1] += Zt[i,pos-j-1]
    # Calcul de P0
    seqAll = ''.join([ "%s"%record.seq for record in seqList ])
    pt[:,0] = [1+Nc(seqAll,'A')-pt[0,1:].sum(),
               1+Nc(seqAll,'C')-pt[1,1:].sum(),
               1+Nc(seqAll,'G')-pt[2,1:].sum(),
               1+Nc(seqAll,'T')-pt[3,1:].sum()]
    # normalisation de Pt
    pt[:,1:] = pt[:,1:] / (Zt.sum() + 4)
    # normalisation de P0
    pt[:,0]  = pt[:,0] / sum(pt[:,0])
    return pt

def get_motives(record, W):
    seq = "%s"%record.seq
    return [ seq[i:i+W] for i in range(len(seq)-W+1) ]

def M_step_v2(seqList,W,Zt):
    # Calcul de Pck
    pt = np.ones((4,W+1))
    for c in ['A','C','G','T']:
        for pos in range(1,W+1):
            for i,record in enumerate(seqList):
                for j,motif in enumerate(get_motives(record,W)):
                    if motif[pos-1] == c:
                        pt[code_dict[c],pos] += Zt[i,j]
    # Calcul de P0
    pt[:,0] = [1+Nc(seqList,'A')-pt[0,1:].sum(),
               1+Nc(seqList,'C')-pt[1,1:].sum(),
               1+Nc(seqList,'G')-pt[2,1:].sum(),
               1+Nc(seqList,'T')-pt[3,1:].sum()]
    # normalisation de Pt
    pt[:,1:] = pt[:,1:] / (Zt.sum() + 4)
    # normalisation de P0
    pt[:,0]  = pt[:,0] / sum(pt[:,0])
    return pt

def MEME(seqList,W,epsilon):
    # init p0
    pt = init_p0_random(seqList,W)
    print "p0 : \n", pt,"\n"
    vr_list = []
    # Iteration until change in p(t) < epsilon
    while len(vr_list) < 2 or abs(vr_list[-2] - vr_list[-1]) > epsilon :
        Zt = E_step(seqList,pt,W)
        pt = M_step_v2(seqList,W,Zt)
        vr_list.append(likelihood_log(seqList,Zt,pt,W))
    return pt,Zt

def get_motif(seqList,indice,pos,W):
    seq = "%s"%seqList[indice].seq
    return seq[pos:pos+W]

def run_meme(handle,W=14):
    seqList = readFasta(handle)
    # print init_p0_random(seqList,W)
    pt,Zt =  MEME(seqList,W,0.001)
    motif_list = []
    print "pt : \n",pt
    print '\n'
    print "Zt : \n",Zt,"\n"
    print "motives positions : ", np.argmax(Zt,axis=1),"\n"
    print "motives : \n"
    for indice,pos in enumerate(np.argmax(Zt,axis=1)):
        motif_list.append(get_motif(seqList,indice,pos,W))
        print motif_list[-1]
    return motif_list

def motif_count(motif_list):
    W = len(motif_list[0])
    m = np.zeros((4,W))
    for i in range(W):
        for c in ['A','C','G','T']:
            m[code_dict[c],i] = Nc(motif_list[:,i],c)
    return m
    
'''
handle = open("hw1_hidden_motif.txt")
motif_list = run_meme(handle)
handle.close()
'''
motif_list = np.array([list('ATGAAAGTTGCGGA'),
                       list('ATCAGAATACCTGA'),
                       list('ATGAAGGGTCCTGG'),
                       list('TTAAATACTACTGA'),
                       list('ATGATAAGTCATCA'),
                       list('ATGGAAGGTACACA'),
                       list('GTAAAAACAACGGA'),
                       list('AGCACAATTATTTA'),
                       list('ATAACAATTACTGA'),
                       list('ATGTAAAGTAAGGA')])

print motif_count(motif_list)
