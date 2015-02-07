#!/usr/bin/python 

# -*- coding: utf-8 -*-

import math
import readline
import numpy as np


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
    for i, l in enumerate(alph):
        alphDict[l] = i
    return alphDict

def key2value(alphDict, keys):
    return np.array([ alphDict.get(k) for k in keys])

def loadA(fname):
    fd = open(fname)
    A = []
    seq = ""
    for line in fd:
        if line[0] == '>' :
            if len(seq) != 0:
                A.append(list(seq))
                seq = ""
        else:
            seq = seq + line.strip()

    if len(seq) != 0:
        A.append(list(seq))
        seq = ""

    fd.close()
    return np.array(A)

# exemple:
A = np.array([[1,2,2,3,1,3],[1,0,2,2,1,3],[1,2,0,2,1,3],[3,2,2,3,1,3]])
seq = [1,2,2,1,3,3]
alph = range(4) # [-,A,B,C]

# Premiere fonction
def count(A, l):
    return sum(A==l)

def n(A, alph):
    return np.array([[count(A[:,i],l) for i in range(len(A[0]))] for l in alph])

def w(A, alph):
    M = len(alph)
    q = len(A)
    return (n(A,alph)+1.)/(M+q) 

W = w(A, alph) 
print 'exemple w:\n', W

# Deuxieme fonction
def S(A, alph):
    q = len(A)
    return np.log2(q) + sum(w(A, alph) * np.log2(w(A, alph)))

print 'exemple S:\n', S(A, alph)

def argmax(A, alph):
    return np.argmax(w(A, alph),axis=0)

print "argmax: \n",argmax(A,alph)

def trace(values, pos):
    pass

# Troisieme fonction
def P0(W,B):
    return np.prod([ W[B[i],i] for i in range(len(B))])

print "P0: \n", P0(W,seq)

def f0(W,B):
    return [ sum(W[b])/len(W[b]) for b in B ]

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

# Appliquer a testseq.txt
alph = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]
print "alphabet: \n", alph
print "length q: \n", len(alph)

seq = loadSeq("test_seq.txt")
print "sequence: \n", seq
alphDict = alph2dict(alph)
alph = key2value(alphDict, alph)
seq  = key2value(alphDict, seq)

A = loadA("Dtrain.txt")
(lign,col) = np.shape(A)
A = key2value(alphDict, A.reshape(lign*col))
A = A.reshape(lign,col)

W = w(A, alph)
print "W0('-'): \n",W[-1,0]
print "S0: \n",S(A, seq[0:len(W[0])])
L = l(W,seq)
print "L:\n", L
