# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np

# données
N = 100
#f = lambda x: x**2
f = lambda x: np.sin(2*x)
#f = lambda x: x**3
datax = -2. + np.random.rand(N)*4.
bruit = np.random.normal(0,0.2,N)       # bruit gaussien
datay = f(datax) + bruit

sigmoid   = lambda a: 1./(1 + np.exp(-a))             
sig_prime = lambda a: sigmoid(a)* (1 - sigmoid(a)) 

# APPRENTISSAGE PAR L'ALGORITHME DE RETROPROPAGATION
T   = 50000                             # nombre total des etapes d'apprentissage
eta = 0.01                              # taux d'apprentissage (learning rate)
W1  = np.random.rand(10,2)-0.5          # premiere couche de poids : 1 neurones d'entrée + biais, 10 neurones cachés
W2  = np.random.rand(11)-0.5            # deuxieme couche de poids : 10 neurones caché + biais, 1 sortie

for i in range(T):
    
    # activite du reseau
    p  = np.random.randint(N)           # choix d'un indice aleatoire parmi N
    x  = np.array([1, datax[p]])        # vecteur des entrées
    a  = np.dot( W1, x)                 # activation des neurones caches
    z  = sigmoid(a)                     # activite des neurones caches

    z  = np.insert(z, 0, 1.)            # ajouter le bias aux neurones cachés
    y  = np.dot(W2, z)                  # calculer la activite du neurone de sortie

    # apprentissage
    target = datay[p]

    delta = y - target                  # erreur delta
    W2 = W2 - eta * delta * z           # regle d'apprentissage pour la deuxieme couche

    delta_h   = sig_prime(a) * delta * W2[1:]   # retropropagation d'erreur
    W1 = W1 - eta * np.outer(delta_h, x)        # regle d'apprentissage pour la premiere couche


# REGRESSION : TRANSFORMATION DES ENTREES PAR LE RESEAU
entrees  = np.arange(-2,2,0.1)
sorties  = []

for xi in entrees:
        
    x = [1, xi]                  # entrée avec biais 
    a = np.dot( W1, x)
    z = sigmoid(a)

    z = np.insert(z, 0, 1)       # neurones cachés avec biais 
    y = np.dot(W2, z)
    
    sorties.append( y )
    

# graphiques
plt.figure(1); plt.clf(); # plt.show()

# l' ensemble des donnees + regression
plt.subplot(111)
plt.scatter( datax, datay, c='r', label='donnees' )
plt.plot(entrees,sorties, 'k-', label = 'regression')
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.draw()

plt.show()

