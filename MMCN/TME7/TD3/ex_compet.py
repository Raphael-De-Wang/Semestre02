# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

# génération de données (deux nuages gaussiens)
N = 100                                          
m1x1, m1x2, s1 =  -1,  5., 1. 
m2x1, m2x2, s2 =   5., 1., 1.

c1x1 = np.random.normal(m1x1, s1, N/2)               
c1x2 = np.random.normal(m1x2, s1, N/2)               

c2x1 = np.random.normal(m2x1, s2, N/2)               
c2x2 = np.random.normal(m2x2, s2, N/2)               

x1  = np.hstack([c1x1, c2x1])             # stocker les coordonees x des deux classes 
x2  = np.hstack([c1x2, c2x2])             # stocker les coordonees y des deux classes

# paramètres de simulation
T = 10000
eta = 0.001

# initisalisation du vecteur de poids
# W = np.random.rand(2,2) - 0.5
n1 = np.random.randint(N)
n2 = np.random.randint(N)
W  = np.array([[x1[n1],x2[n1]],[x1[n2],x2[n2]]])

plt.figure(1); plt.clf(); plt.show()   # initialisation graphique

for i in range(T):
    
    p = np.random.randint(N)
    x = np.array([x1[p], x2[p]])
    y = np.dot(W, x)

    # ajouter l'apprentissage compétitif pour la matrice de W
    # afin de trouver l'index de la sortie maximale, utiliser np.argmax(y)
    i_max = np.argmax(y)
    W[i_max] += eta*(x - W[i_max]) 

    if i%100==0:

        print('Iteration %d' % i)
        plt.subplot(111); plt.cla()
        plt.scatter( x1, x2 )  
        plt.plot([0, W[0,0]], [0, W[0,1]], '-r', lw=3)
        plt.plot([0, W[1,0]], [0, W[1,1]], '-g', lw=3)
        plt.axhline(0, ls=':', color='k')                    
        plt.axvline(0, ls=':', color='k')                    
        plt.xlabel("x1")
        plt.ylabel("x2")
        plt.axis('scaled')                                    
        plt.draw()

        #sleep(0.5)
