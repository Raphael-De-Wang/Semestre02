#!/usr/bin/python 
# -*- coding: utf-8 -*-

import math
import readline
import numpy as np
import matplotlib.pyplot as plt

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

D = np.array([[0,1,1,2,0,2],
              [0,3,1,1,0,2],
              [0,1,3,1,0,2],
              [2,1,1,2,0,2]])

def count(D,a,b,i,j):
    return sum(( D[:,i] == a ) * (D[:,j] == b))

def Nij(D,a,b):
    length = len(D[0])
    N = np.zeros((length,length))
    for i in range(length):
        for j in range(length):
            N[i,j] = count(D,a,b,i,j)
    return N

N = Nij(D,1,2)
print Nij(D,1,0)+Nij(D,1,1)+Nij(D,1,2)+Nij(D,1,3)
print Nij(D,0,2)+Nij(D,1,2)+Nij(D,2,2)+Nij(D,3,2)

def Ni(D,l):
    return np.array([sum(D[:,i]==l) for i in range(len(D[0]))])

print Ni(D,1)
print Ni(D,2)

def w(N,M,q):
    return (N + (1./q))/(M+q)
    
print w(N,4,4)

def m(D,alph,M,q):
    length = len(D[0])
    Mij    = np.zeros((length,length))
    for a in alph:
        for b in alph:
            N = Nij(D,a,b)
            Na= Ni(D,a)
            Nb= Ni(D,b)
            W = w(N,M,q)
            Wa= w(Na,M,q)
            Wb= w(Nb,M,q)
            Mij = Mij + W*(np.log2(W/np.outer(Wa,Wb)))
    return Mij
    
print m(D,[0,1,2,3],4,4)

alph = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-"]
alphDict,valDict = alph2dict(alph)
alph = key2value(alphDict, alph)

D = loadD("TD1/Dtrain.txt")
(lign,col) = np.shape(D)
D = key2value(alphDict, D.reshape(lign*col))
D = D.reshape(lign,col)

M = len(D)
q = len(alph)

a = 0
b = 1
Wa = w(Ni(D,a),M,q)
print "Wa: \n",Wa
np.savetxt("Wa.txt",Wa,fmt="%.4f")
N = Nij(D,a,b)
print "N: \n",N
np.savetxt("N.txt",N,fmt="%.4f")
Wab = w(N,M,q)
print "Wab: \n",Wab
np.savetxt("Wab.txt",Wab,fmt="%.4f")
'''
Mij = m(D,alph,M,q)
print Mij
'''
