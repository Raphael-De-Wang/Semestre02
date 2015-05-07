# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np
from donnees import *

# CREATION DE L'ENSEMBLE DES DONNEES 
#
# créer un ensemble de données (2 classes en 2D)

N=200

datax1,datax2,labels = gen_2class(N)

sigmoid   = lambda a: 1./(1 + np.exp(-a))   # definition d'une function sigmoide

# initialiser le vecteur des poids 
w = np.random.rand(2) - 0.5

# CLASSIFICAITON (sans apprentissage)
classe = []                                 # initialisation du tableau pour garder les resultats de classification

for i in range(N):
    
    motif = np.array([datax1[i], datax2[i]])    # i-ème motif 
    a     = np.dot( w, motif )                  # activation du neurone de sortie
    y     = sigmoid(a)                          # activite du neurone de sortie

    # classification
    if y>0.5: 
        classe.append('r')
    else:
        classe.append('b')


# figures
plt.figure(1); plt.clf(); # plt.show()                

# l'ensemble des donnees
plt.subplot(221)
plt.scatter( datax1, datax2, c = labels )            # donnees
plt.axhline(0, ls=':', color='k')                    # ligne horizontale pointillee (':') noire ('k')
plt.axvline(0, ls=':', color='k')                    # ligne verticale pointillee (':') noire ('k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')                                    
plt.title("L'ensemble des donnees")
plt.draw()

# resultat de classification (sans apprentissage)
plt.subplot(222)
plt.scatter( datax1, datax2, c = classe )
plt.quiver( 0, 0, w[0]/np.linalg.norm(w), w[1]/np.linalg.norm(w), angles='xy',scale_units='xy',scale=1 ) # vecteur normal
ax1, bx1 = -2, 2
ax2, bx2 =  - w[0]*ax1 / w[1], - w[0]*bx1 / w[1]
plt.plot([ax1, bx1], [ ax2, bx2], 'g' )                 # hyperplan separateur
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')
plt.title('classification')
plt.draw()

# fonction d'activation (sigmoide)
plt.subplot(223)
a = np.arange(-10,10,0.1)
plt.plot( a, sigmoid(a) )
plt.axhline(0.5, ls='--', color='r')
plt.xlabel("a")
plt.ylabel("sigmoid(a)")
plt.title ("fonction d'activation g(a)")
plt.draw()

# plt.show()
