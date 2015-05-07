# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np

# CREATION DE L'ENSEMBLE DES DONNEES (3 classes avec la distribution Gaussienne)
N = 150
m1x, m1y, s1 = 1, 10, 1
m2x, m2y, s2 = 5, 5, 1
m3x, m3y, s3 = 1, 1, 1

c1x1 = np.random.normal(m1x, s1, N/3)               # coordonnees x des points dans la classe 1
c1x2 = np.random.normal(m1y, s1, N/3)               # coordonnees y des points

c2x1 = np.random.normal(m2x, s2, N/3)               # coordonnees x des points dans la classe 2
c2x2 = np.random.normal(m2y, s2, N/3)               # coordonnees y des points

c3x1 = np.random.normal(m3x, s3, N/3)               # coordonnees x des points dans la classe 3
c3x2 = np.random.normal(m3y, s3, N/3)               # coordonnees y des points

datax1  = np.hstack([c1x1, c2x1, c3x1])             # stocker les coordonees x des deux classes 
datax2  = np.hstack([c1x2, c2x2, c3x2])             # stocker les coordonees y des deux classes

labels = ['r']*(N/3) + ['g']*(N/3) + ['b']*(N/3)    # etiquettes de motifs 
                                                    # 'r' - classe 1, 'g' - classe 2, 'b' - classe 3

index  = np.random.permutation(range(N))            # permuter aleatoirement les indices
datax1  = datax1[index]                             # permuter les x    
datax2  = datax2[index]                             # permuter les y
labels  = [labels[i] for i in index]                # permuter les etiquettes 
                                                    # (traitement special car 'labels' est un tableau des objets)

sigmoid   = lambda a: 1./(1 + np.exp(-a))           # definition d'une function sigmoide
sig_prime = lambda a: sigmoid(a)* (1 - sigmoid(a))  # derivee de la fonction sigmoide

# APPRENTISSAGE PAR LA DESCENTE DU GRADIENT
T = 1000                                        # nombre d'itérations
eta = 0.1                                       # taux d'apprentissage (learning rate)
W = np.random.rand(3,3)-0.5                     # initialisation aleatoire des poids synaptiques

for i in range(T):
    
    p = np.random.randint(N)                    # choix d'un indice aleatoire parmi N
    x = np.array([1., datax1[p], datax2[p]])    # motif  x(p)
    a = np.dot( W, x)                           # activation du neurone de sortie
    y = sigmoid(a)                              # activite du neurone de sortie

    # apprentissage
    if labels[p] == 'r':
        target = np.array([1,0,0])
    elif labels[p] == 'g':
        target = np.array([0,1,0])
    else:
        target = np.array([0,0,1])
 
    delta = sig_prime(a) * ( y - target )
    W = W - eta * np.outer(delta, x)

# CLASSIFICATION PAR l'APPRENTISSAGE PAR LA DESCENTE DU GRADIENT
classe = []                                 # initialisation du tableu pour garder les resultats de classification
for i in range(N):
    
    # activite du reseau
    x = [ 1, datax1[i], datax2[i] ]         # i-ème motif 
    a = np.dot( W, x )                      # activation du neurone de sortie
    y = sigmoid(a)                          # activite du neurone de sortie

    # classification
    ind_max = np.argmax(y)

    if ind_max == 0: 
        classe.append('r')
    elif ind_max == 1:
        classe.append('g')
    else:
        classe.append('b')


# graphiques
plt.figure(1); plt.clf(); # plt.show()                

# l'ensemble des donnees
plt.subplot(211)
plt.scatter( datax1, datax2, c = labels )            # donnees
plt.axhline(0, ls=':', color='k')                    # ligne horizontale pointillee (':') noire ('k')
plt.axvline(0, ls=':', color='k')                    # ligne verticale pointillee (':') noire ('k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')                                    
plt.title("L'ensemble des donnees")
plt.draw()

# resultat de classification 
plt.subplot(212)
plt.scatter( datax1, datax2, c = classe )
ax1, bx1 = -2., 6.
ax2rg = -(W[0,1]-W[1,1])/(W[0,2]-W[1,2])*ax1 - (W[0,0]-W[1,0])/(W[0,2]-W[1,2])
bx2rg = -(W[0,1]-W[1,1])/(W[0,2]-W[1,2])*bx1 - (W[0,0]-W[1,0])/(W[0,2]-W[1,2])
ax2rb = -(W[0,1]-W[2,1])/(W[0,2]-W[2,2])*ax1 - (W[0,0]-W[2,0])/(W[0,2]-W[2,2])
bx2rb = -(W[0,1]-W[2,1])/(W[0,2]-W[2,2])*bx1 - (W[0,0]-W[2,0])/(W[0,2]-W[2,2])
ax2gb = -(W[1,1]-W[2,1])/(W[1,2]-W[2,2])*ax1 - (W[1,0]-W[2,0])/(W[1,2]-W[2,2])
bx2gb = -(W[1,1]-W[2,1])/(W[1,2]-W[2,2])*bx1 - (W[1,0]-W[2,0])/(W[1,2]-W[2,2])
plt.plot([ax1, bx1], [ ax2rg, bx2rg], 'k' )                 # hyperplan separateur
plt.plot([ax1, bx1], [ ax2rb, bx2rb], 'k' )                 # hyperplan separateur
plt.plot([ax1, bx1], [ ax2gb, bx2gb], 'k' )                 # hyperplan separateur
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.axis('equal')
plt.xlim(-4,8)
plt.ylim(-2,14)
plt.xlabel("x1")
plt.ylabel("x2")
plt.title('classification')
plt.draw()

plt.show()                
