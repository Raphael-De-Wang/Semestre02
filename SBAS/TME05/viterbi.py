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

def persistHMM(fname,data):
    np.save(fname,data)

def loadHMM(fname):
    return np.load(fname)
    
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

    return [etatsCaches,observs]

# Viterbi
def viterbi(observs,PI0,P,E):
    v = np.array([PI0])*E[:,observs[0]]
    for o in observs[1:]:
        v = np.vstack((v,max(v[-1])*P[np.argmax(v[-1])]*E[:,observs[o]]))
    return (max(v[-1]),np.argmax(v,axis=1))

# Estimation des parametres :
def sampling(PI0,P,E,T,nombre):
    return np.array([ simulation(PI0,P,E,T) for i in xrange(nombre) ])

def estiP(etatsCaches,N):
    P = np.ones((N,N)) # au lieu de np.zeros() pour evider " RuntimeWarning: invalid value encountered in divide "
    for chain in etatsCaches:
        for i in range(len(chain)-1):
            P[chain[i],chain[i+1]] = P[chain[i],chain[i+1]] + 1
    return np.transpose(P.T/np.sum(P,axis=1))
    
def estiE(etatsCaches,observs,N,M):
    E = np.ones((N,M)) # au lieu de np.zeros() pour evider " RuntimeWarning: invalid value encountered in divide "
    for i in range(len(etatsCaches)):
        for j in range(len(etatsCaches[0])):
            E[etatsCaches[i,j],observs[i,j]] = E[etatsCaches[i,j],observs[i,j]] + 1
    return np.transpose(E.T/np.sum(E,axis=1))

#  Viterbi training :
def viterbiTraining(observs,PI0,N,M,nomIter):
    # initialiser un point aleatoire dans le convergence. Pour simplifier, je prend uniform
    P = np.ones((N,N))/N
    E = np.ones((N,M))/M
    for i in xrange(nomIter):
        etatsCaches = None
        for obs in observs:
            prob,X = viterbi(obs,PI0,P,E)
            if etatsCaches == None:
               etatsCaches = np.array([X])
            else:
                etatsCaches = np.append(etatsCaches,[X],axis=0)
        E = estiE(etatsCaches,observs,N,M)
        P = estiP(etatsCaches,N)
    return (P,E)

# 3c) (Viterbi training deuxime version)
def viterbiTraining2(observs,PI0,N,M,nomIter):
    pass

def logVraisemblance():
    pass
    
# Q1
t = 200
'''
Xt,Yt=simulation(pi0,p,e,t)
print Xt,Yt
'''

# Q2
# print viterbi(Xt,pi0,p,e)

# Q3
# print sampling(pi0,p,e,t,10)

# Q3.A
'''
print estiP(np.array([[0,0,1],[0,1,0]]),2)
print estiE(np.array([[0,0,1],[0,1,0]]),np.array([[0,0,1],[1,1,2]]),2,3)
'''

# Q3.B
'''
obs = sampling(pi0,p,e,t,100)[:,1]
print viterbiTraining(obs,pi0,2,2,20)
'''

# Q3.C