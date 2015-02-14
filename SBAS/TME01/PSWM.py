#!/usr/bin/python 

# -*- coding: utf-8 -*-

import math
import readline
import numpy as np
import matplotlib.pyplot as plt

# np.set_printoptions(precision=2, linewidth=320)

# example pour tester
def loadSeq(fname):
    fd = open(fname)
    seq = ""
    for line in fd:
        if line[0] == '>' :
            continue
        else:
            seq = seq + line.strip()
    fd.close()
    return list(seq)

def alph2dict(alph):
    alphDict = {}
    valDict  = {}
    for i, l in enumerate(alph):
        alphDict[l] = i
        valDict["%d"%i] = l
    return alphDict,valDict

def key2value(alphDict, keys):
    return np.array([ alphDict.get(k) for k in keys])

def loadD(fname):
    fd = open(fname)
    D = []
    seq = ""
    for line in fd:
        if line[0] == '>' :
            if len(seq) != 0:
                D.append(list(seq))
                seq = ""
        else:
            seq = seq + line.strip()

    if len(seq) != 0:
        D.append(list(seq))
        seq = ""

    fd.close()
    return np.array(D)

# exemple:
D = np.array([[0,1,1,2,0,2],
              [0,3,1,1,0,2],
              [0,1,3,1,0,2],
              [2,1,1,2,0,2]])
seq = [0,1,1,0,0,2]
alph = range(4) # [A,B,C,-]

# Premiere fonction
def count(D, l):
    return sum(D==l)

def n(D, alph):
    return np.array([[count(D[:,i],l) for i in range(len(D[0]))] for l in alph])
    
print 'exemple n:\n', n(D, alph)

def w(D, alph):
    M = len(D)
    q = len(alph)
    return (n(D,alph)+1.)/(M+q) 

W = w(D, alph) 
print 'exemple w:\n', W

# Deuxieme fonction
def S(q, W):
    return np.log2(q) + sum(W * np.log2(W))

print 'exemple S:\n', S(len(alph), W)

def argmax(L):
    return np.argmax(L,axis=0)

print "argmax: \n",argmax(D)

def traceS(valeurs, tname):
    fig = plt.figure()
    plt.xlabel("Position")
    plt.ylabel("Entropy")
    plt.plot(valeurs)
    plt.savefig(tname)

# Troisieme fonction
def P0(W,B):
    return np.prod([ W[B[i],i] for i in range(len(B))])

print "P0: \n", P0(W,seq)

def f0(W,B):
    return np.array([ sum(W[b])/len(W[b]) for b in B ])

print "f0:\n",f0(W,alph)

# Quatrieme fonction
def l(W,B):
    L = []
    for j in range(len(B)-len(W[0])+1):
        b = B[j:j+len(W[0])]
        F0 = f0(W,b)
        L.append(sum(np.log2([W[b[i],i]/F0[i] for i in range(len(F0))])))
    return np.array(L)

print "l: \n", l(W,seq)
print "==== test end ====\n"

def traceL(valeurs, tname):
    fig = plt.figure()
    plt.xlabel("Position")
    plt.ylabel("Log-Vraisemblance")
    plt.plot(valeurs)
    plt.savefig(tname)

# Appliquer a testseq.txt
alph = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]

seq = loadSeq("TD1/test_seq.txt")
alphDict,valDict = alph2dict(alph)
alph = key2value(alphDict, alph)
seq  = key2value(alphDict, seq)

D = loadD("TD1/Dtrain.txt")
(lign,col) = np.shape(D)
print "D shape:", lign, col
D = key2value(alphDict, D.reshape(lign*col))
D = D.reshape(lign,col)

N = n(D, alph)
print "N \n", N
print "N Shape:\n", np.shape(N)
np.savetxt("N.txt", N,fmt="%.2f")

W = w(D, alph)
print "W: \n",W
print "W Shape: \n",np.shape(W)
np.savetxt("W.txt", W,fmt="%.2f")

entropy = S(len(alph), W)
print "S0: \n",entropy
print "S0 Length: \n",np.shape(entropy)
np.savetxt("S0.txt",entropy,fmt="%.2f")

modele_nul = f0(W,alph)
print "f0:\n",modele_nul
print "f0 Length:\n",np.shape(modele_nul)
np.savetxt("f0.txt",modele_nul,fmt="%.2f")

L = l(W,seq)
print "L:\n", L
print "L Length:\n", np.shape(L)
np.savetxt("L.txt",L,fmt="%.2f")

pMax = argmax(L)
print "pMax: \n", pMax
print "pMax sequence: \n", key2value(valDict,np.array([ "%d"%i for i in seq[pMax:pMax+col]]))
print "Max Log Vraisemblance: ", l(W,seq[pMax:pMax+col])

traceS(entropy,"traceS")
traceL(L, "traceL")

