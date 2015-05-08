# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# CREATION DE L'ENSEMBLE DES DONNEES (2 classes)
N = 150
m1x, m1y, s1 = 10, 10, 1
m2x, m2y, s2 = 5, 5, 1
m3x, m3y, s3 = 1, 1, 1

c1x1 = np.random.normal(m1x, s1, N/3)
c1x2 = np.random.normal(m1y, s1, N/3)

c2x1 = np.random.normal(m2x, s2, N/3)
c2x2 = np.random.normal(m2y, s2, N/3)

c3x1 = np.random.normal(m3x, s3, N/3)
c3x2 = np.random.normal(m3y, s3, N/3)

datax1  = np.hstack([c1x1, c2x1, c3x1])
datax2  = np.hstack([c1x2, c2x2, c3x2])

labels = ['r']*(N/3) + ['b']*(N/3) + ['r']*(N/3)    # etiquettes de points dans l'ensemble des donnees
                                                    # 'r' - classe 1, 'b' - classe 2

index   = np.random.permutation(range(N))   # permuter aleatoirement les indices
datax1  = datax1[index]
datax2  = datax2[index]
labels  = [labels[i] for i in index]

sigmoid   = lambda a: 1./(1 + np.exp(-a))             
sig_prime = lambda a: sigmoid(a)* (1 - sigmoid(a)) 

# APPRENTISSAGE PAR L'ALGORITHME DE RETROPROPAGATION
T   = 10000                             # nombre total des etapes d'apprentissage
eta = 0.1                               # taux d'apprentissage (learning rate)
D   = 2                                 # dimension de l'espace des entrees
H   = 10                                # nombre de neurones caches
W1  = np.random.rand(10,3) - 0.5        # premiere couche de poids
W2  = np.random.rand(11) - 0.5          # deuxieme couche de poids

for i in range(T):
    
    # activite du reseau
    # activite du reseau
    p  = np.random.randint(N)                   # choix d'un indice aleatoire parmi N
    x  = np.array([1, datax1[p], datax2[p]])    # vecteur des entrées
    a  = np.dot( W1, x)                         # activation des neurones caches
    z  = sigmoid(a)                             # activite des neurones caches

    z  = np.insert(z, 0, 1.)                    # ajouter le bias aux neurones cachés
    y  = np.dot(W2, z)                          # calculer la activite du neurone de sortie

    # apprentissage
    if labels[p] == 'r':
        target = 1
    else: 
        target = 0
    
    delta = y - target                          # erreur delta
    W2 = W2 - eta * delta * z                   # regle d'apprentissage pour la deuxieme couche

    delta_h   = sig_prime(a) * delta * W2[1:]   # retropropagation d'erreur
    W1 = W1 - eta * np.outer(delta_h, x)        # regle d'apprentissage pour la premiere couche
       
# CLASSIFICATION 
classe = []
for i in range(N):
    
    # activite du reseau
    x = np.array([1, datax1[i], datax2[i]])     # i-ème motif
    a = np.dot( W1, x)                          # activation des neurones caches
    z = sigmoid(a)                              # activite des neurones caches

    z = np.insert(z, 0, 1.)                     # ajouter le bias
    y = np.dot(W2, z)                           # calculer la activite du neurone de sortie

    # classification
    if y>0.5: 
        classe.append('r')
    else:
        classe.append('b')

# graphiques
plt.figure(1); plt.clf(); plt.show()                

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

# resultat de classification 
plt.subplot(222)
plt.scatter( datax1, datax2, c = classe )
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')
plt.title('classification')
plt.draw()

# TRANSFORMATION D'ENTREES PAR LE RESEAU
x1range  = np.arange(-5,15,0.1)
x2range  = np.arange(-5,15,0.1)
y = np.zeros( (len(x2range), len(x1range)) )
for i, x2 in enumerate(x2range):
    for j, x1 in enumerate(x1range):
        
        x = np.array([1, x1, x2])
        a = np.dot( W1, x)
        z = sigmoid(a)

        z = np.insert(z, 0, 1)       # bias 
        y[i,j] = np.dot(W2, z)

# figure : transformation d'entrees
X1, X2 = np.meshgrid(x1range, x2range)
ax = plt.subplot(223, projection='3d')
ax.plot_surface(X1, X2, y, rstride=5, cstride=5, cmap='jet')
plt.xticks([-5,15])
plt.yticks([-5,15])
plt.xlabel("x1")
plt.ylabel("x2")
ax.set_zlabel("y")
plt.title("Transformation des entrees")
plt.draw()

plt.subplot(224)
plt.contourf(X1, X2, y, 50)
plt.xlabel("x1")
plt.ylabel("x2")
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.axis('scaled')
plt.title("Transformation des entrees")
plt.draw()

plt.show()
