#!/bin/python

import numpy as np


# donnee initiale
F="FAIR"
U="UNFAIR"
H="HEAD"
T="TAIL"

S=np.array([F,U])
O=np.array([H,T])
p=np.array([[0.99,0.01],[0.05,0.95]])
e=np.array([[0.5,0.5],[0.9,0.1]])
pi0=np.array([1.,0])

# Function basique
def accumulate(M):
    if len(np.shape(M)) == 1:
        return np.cumsum(M)
    return np.cumsum(M,axis=1)

def tirageAleatoire(outPutSize=None):
    return np.random.uniform(size=outPutSize)

def translate(valueList,labels):
    return np.array([ labels[v] for v in valueList ])

def choisiEtat(probDistrb,prob):
    for index,seuil in enumerate(probDistrb):
        if seuil > prob:
            return index

# Simulation :
def simulation(PI0,P,E,T):
    aPI0 = accumulate(PI0)
    aP   = accumulate(P)
    aE   = accumulate(E)

    etatsCaches = np.array([])
    observs     = np.array([])

    etatsCaches = np.append(etatsCaches,choisiEtat(aPI0,tirageAleatoire()))
    observs     = np.append(observs,choisiEtat(aE[etatsCaches[0]],tirageAleatoire()))

    for i in xrange(T-1):
        etatsCaches = np.append(etatsCaches,choisiEtat(aP[etatsCaches[-1]],tirageAleatoire()))
        observs     = np.append(observs,choisiEtat(aE[etatsCaches[-1]],tirageAleatoire()))

    return (etatsCaches,observs)

# Viterbi
def viterbi(observs,PI0,P,E):
    v = np.array([PI0])*E[:,observs[0]]
    for o in observs[1:]:
        v = np.vstack((v,max(v[-1])*P[np.argmax(v[-1])]*E[:,observs[o]]))
    return (max(v[-1]),np.argmax(v,axis=1))

t = 10
Xt,Yt=simulation(pi0,p,e,t)
print Xt,Yt
print viterbi(Xt,pi0,p,e)
