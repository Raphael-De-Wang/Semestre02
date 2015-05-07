# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np


def gen_2class(N=150): # la taille de l'ensemble
    # CREATION DE L'ENSEMBLE DES DONNEES (3 classes avec la distribution Gaussienne)
    m1x, m1y, s1 = 1., 10., 1.                          # moyenne (x,y) et largeur (s) de la classe 1
    m2x, m2y, s2 = 5., 5., 1.                           # moyenne et largeur de la classe 2
    
    c1x1 = np.random.normal(m1x, s1, N/2)               # coordonnees x des points dans la classe 1
    c1x2 = np.random.normal(m1y, s1, N/2)               # coordonnees y des points
    
    c2x1 = np.random.normal(m2x, s2, N/2)               # coordonnees x des points dans la classe 2
    c2x2 = np.random.normal(m2y, s2, N/2)               # coordonnees y des points
    
    datax1  = np.hstack([c1x1, c2x1])             # stocker les coordonees x des deux classes 
    datax2  = np.hstack([c1x2, c2x2])             # stocker les coordonees y des deux classes
    
    labels = ['r']*(N/2) + ['g']*(N/2)                  # etiquettes de points dans l'ensemble des donnees
    # 'r' - classe 1, 'g' - classe 2

    index  = np.random.permutation(range(N))            # permuter aleatoirement les indices
    datax1  = datax1[index]                             # permuter les x    
    datax2  = datax2[index]                             # permuter les y
    labels  = [labels[i] for i in index]                # permuter les etiquettes 
    # (traitement special car 'labels' est un tableau des objets)
    return datax1,datax2,labels

def gen_3class(N=150): # la taille de l'ensemble
    # CREATION DE L'ENSEMBLE DES DONNEES (3 classes avec la distribution Gaussienne)
    m1x, m1y, s1 = 1., 10., 1.                          # moyenne (x,y) et largeur (s) de la classe 1
    m2x, m2y, s2 = 5., 5., 1.                           # moyenne et largeur de la classe 2
    m3x, m3y, s3 = 1., 1., 1.                           # moyenne et largeur de la classe 3
    
    c1x1 = np.random.normal(m1x, s1, N/3)               # coordonnees x des points dans la classe 1
    c1x2 = np.random.normal(m1y, s1, N/3)               # coordonnees y des points
    
    c2x1 = np.random.normal(m2x, s2, N/3)               # coordonnees x des points dans la classe 2
    c2x2 = np.random.normal(m2y, s2, N/3)               # coordonnees y des points
    
    c3x1 = np.random.normal(m3x, s3, N/3)               # coordonnees x des points dans la classe 3
    c3x2 = np.random.normal(m3y, s3, N/3)               # coordonnees y des points
    
    datax1  = np.hstack([c1x1, c2x1, c3x1])             # stocker les coordonees x des deux classes 
    datax2  = np.hstack([c1x2, c2x2, c3x2])             # stocker les coordonees y des deux classes
    
    labels = ['r']*(N/3) + ['g']*(N/3) + ['b']*(N/3)    # etiquettes de points dans l'ensemble des donnees
    # 'r' - classe 1, 'g' - classe 2, 'b' - classe 3

    index  = np.random.permutation(range(N))            # permuter aleatoirement les indices
    datax1  = datax1[index]                             # permuter les x    
    datax2  = datax2[index]                             # permuter les y
    labels  = [labels[i] for i in index]                # permuter les etiquettes 
    # (traitement special car 'labels' est un tableau des objets)
    return datax1,datax2,labels

def plot_data(datax1,datax2):
    # graphiques
    plt.figure(1); plt.clf(); 
    
    # l'ensemble des donnees
    plt.subplot(111)
    plt.scatter( datax1, datax2, c = labels )            # donnees
    plt.axhline(0, ls=':', color='k')                    # ligne horizontale pointillee (':') noire ('k')
    plt.axvline(0, ls=':', color='k')                    # ligne verticale pointillee (':') noire ('k')
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.axis('scaled')                                    
    plt.title("L'ensemble des donnees")
    plt.draw()
    
    plt.show()                
