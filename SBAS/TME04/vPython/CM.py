#!/bin/python

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

#              TO SL,   Basal, Luminal  FROM
P159 = np.array([[0.58, 0.35,  0.07], # SL
                 [0.01, 0.99,  0.00], # Basal
                 [0.04, 0.49,  0.47]])# Luminal

#              TO SL,   Basal, Luminal  FROM
P149 = np.array([[0.61, 0.09,  0.30], # SL
                 [0.01, 0.90,  0.08], # Basal
                 [0.01, 0.00,  0.99]])# Luminal

P149[1] = P149[1]/sum(P149[1])

# 1) Traduire les graphes de transition en matrices de transition, P 149 et P 159 . Creer les matrice en R. Determiner les etats stationnaires PI149 et PI159 (commandes matrix, eigen). Est-ce que les valeurs que vous trouvez sont consistentes avec les valeurs experimentales en Fig. 1?
# Attention: Les probabilites de transition indiquees dans le graphe de SUM149 ne sont pas normalisees. Il faut corriger cette petite imprecision dans votre implementation.

def eigen(P):
    (valeurProp, vecteurProp) = linalg.eig(P.T)
    pi = abs(vecteurProp[:,np.argmax(abs(valeurProp))])
    return pi / sum(pi)

P159s = eigen(P159)
P149s = eigen(P149)
'''
print P149s
'''

#                 SL,    Basal, Luminal
PI0 = np.array([0.019, 0.973, 0.0062])
PI0 = np.array([0.039, 0.033, 0.928])

# 2) Determiner les matrices de transition P 149,159 pour n iterations avec n = 2, 4, 8, 16, 32, 64. Comparer les matrices avec Pi 149,159 .
def iter_approcher(PI0,P,PI,eps=0.001):
    Pn = PI0.dot(P)
    distance = np.inf
    i = 0
    while np.linalg.norm(PI-Pn) > eps:
        i = i + 1
        Pn = Pn.dot(P)
    return (i, Pn)

'''
(i,Pn) = iter_approcher(PI0,P159,P159s)
print Pn
print i
print P159s
'''

'''
(i,Pn) = iter_approcher(PI0,P159,P159s)
print Pn
print i
print P159s
'''

def vecPuis(v,m,n):
    for i in xrange(n):
        v = np.dot(v,m)
    return v

'''
for n in range(1,7):
    print vecPuis(PI0, P159, 2**n)
'''

# 3) Ecrire une fonction pour determiner la probabilite ( PI(t) ) t=0:T d'une chaine de Markov a trois etats (trajectoire moyenne).
# Appliquer pour P 149 et P 159 , pour les conditions initiales PI(0): {(0.998, 0.001, 0.001), (0.001, 0.998, 0.001), (0.001, 0.001, 0.998)} et pour T = 20 iterations de la chaine de Markov.
# Visualiser graphiquement les resultats.

PI0 = np.array([(0.998, 0.001, 0.001), (0.001, 0.998, 0.001), (0.001, 0.001, 0.998)])

def expand(PI0,P,T):
    pi0 = np.array([PI0])
    for t in xrange(T):
        pi0 = np.vstack((pi0,[pi0[-1].dot(P)]))
    return pi0

T = 20
# print expand(PI0, P159, T)
# print expand(PI0, P149, T)

def curveSBL(S,B,L):
    fig = plt.figure()
    plt.plot(S)
    plt.plot(B)
    plt.plot(L)
    plt.show()

def etat_suivant(P,etat_actual):
    accP = np.add.accumulate(P,1)
    alea = np.random.uniform()
    for i in range(len(accP[etat_actual])):
        if alea <= accP[etat_actual,i]:
            return i
    
def CM(etat_inital,T,P):
    chaine = [etat_inital]
    for t in xrange(T):
        chaine.append(etat_suivant(P,chaine[-1]))
    return chaine
        
def populationCM(etat_inital,N,T,P):
    return np.array([ CM(etat_inital,T,P) for n in xrange(N) ])

N       = 1000
SL      = 0
Basal   = 1
Luminal = 2
P       = P159
# P     = P149
chaineSL = populationCM(SL, N, T, P)
chaineBa = populationCM(Basal, N, T, P)
chaineLu = populationCM(Luminal, N, T, P)

def rate(popChains, etat):
    popul = len(popChains)
    return sum(popChains == etat, 0 ) * 1. / popul
'''
rateSL   = rate(chaineSL, SL) # sum(chaineSL == SL, 0 ) * 1. / N
rateBa   = rate(chaineSL, Basal) # sum(chaineSL == Basal, 0 ) * 1. / N
rateLu   = rate(chaineSL, Luminal) # sum(chaineSL == Luminal, 0 ) * 1. / N

curveSBL(rateSL,rateBa,rateLu)
'''

def prob_traj(chaine,P):
    prob = 0.
    for i in range(len(chaine)-1):
        # print P[chaine[i],chaine[i+1]]
        prob = prob + np.log(P[chaine[i],chaine[i+1]])
    return prob
    
# print prob_traj(chaineSL[0],P)

probsSL = [ prob_traj(s,P) for s in chaineSL ]
probsBa = [ prob_traj(s,P) for s in chaineBa ]
probsLu = [ prob_traj(s,P) for s in chaineLu ]

# (y,x) = np.histogram(probsSL)
fig = plt.figure()

axSL = fig.add_subplot(1, 1, 1)
axBa = fig.add_subplot(1, 1, 1)
axLu = fig.add_subplot(1, 1, 1)
bins = 25

axSL.hist(probsSL,bins)
axBa.hist(probsBa,bins)
axLu.hist(probsLu,bins)

plt.show()

