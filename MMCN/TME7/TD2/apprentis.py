# -*-coding:Latin-1 -*
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# CREATION DE L'ENSEMBLE DES DONNEES (2 classes avec la distribution Gaussienne)
N = 100                                     # la taille de l'ensemble
m1x, m1y, s1 =   5.,  5., 2.                # moyenne (x,y) et largeur (s) de la classe 1
m2x, m2y, s2 =   1.,  1., 2.                # moyenne et largeur de la classe 2

c1x1 = np.random.normal(m1x, s1, N/2)       # coordonnees x des points dans la classe 1
c1x2 = np.random.normal(m1y, s1, N/2)       # coordonnees y des points
    
c2x1 = np.random.normal(m2x, s2, N/2)       # coordonnees x des points dans la classe 2
c2x2 = np.random.normal(m2y, s2, N/2)       # coordonnees y des points

datax1  = np.hstack([c1x1, c2x1])           # stocker les coordonees x des deux classes
datax2  = np.hstack([c1x2, c2x2])           # stocker les coordonees y des deux classes

labels = ['r']*(N/2) + ['b']*(N/2)          # etiquettes de points dans l'ensemble des donnees
                                            # 'r' (rouge) - classe 1, 'b' (bleu) - classe 2

index  = np.random.permutation(range(N))    # permuter aleatoirement les indices 
datax1  = datax1[index]                     # permuter les x    
datax2  = datax2[index]                     # permuter les y
labels  = [labels[i] for i in index]        # permuter les etiquettes 
                                            # (traitement special car 'labels' est un tableau des objets)

sigmoid   = lambda a: 1./(1 + np.exp(-a))           # definition d'une function sigmoide
sig_prime = lambda a: sigmoid(a)* (1 - sigmoid(a))  # derivee de la fonction sigmoide

print "# APPRENTISSAGE PAR LA DESCENTE DU GRADIENT"
T = 1000                                    # nombre d'itérations
eta = 0.1                                   # taux d'apprentissage (learning rate)
w = np.random.rand(3)-0.5                   # initialisation aleatoire des poids synaptiques

for i in range(T):
    
    p = np.random.randint(N)                # choix d'un indice aleatoire parmi N
    x = np.array([1, datax1[p], datax2[p]])    # motif  x(p)
    a = np.dot( w, x)                       # activation du neurone de sortie
    y = sigmoid(a)                          # activite du neurone de sortie

    # conversion des étiquettes à 0 et 1
    if labels[p] == 'r':
        target = 1.
    else: 
        target = 0.
        
    # delta rule
    delta = sig_prime(a) * ( y - target )   # erreur delta
    w = w - eta * delta * x                 # descente du gradient

      
print "# CLASSIFICATION"
classe = []                      # initialisation du tableau pour garder les resultats de classification
for i in range(N):
    
    x = [ 1, datax1[i], datax2[i] ]            # i-ème motif 
    a = np.dot( w, x )                      # activation du neurone de sortie
    y = sigmoid(a)                          # activite du neurone de sortie

    # classification
    if y>0.5: 
        classe.append('r')
    else:
        classe.append('b')


print "# graphiques"
plt.figure(1); plt.clf(); # plt.show()                

# l'ensemble des donnees
plt.subplot(221)
plt.scatter( datax1, datax2, c = labels )            # donnees
plt.axhline(0, ls=':', color='k')                    # ligne horizontale pointillee (':') noire ('k')
plt.axvline(0, ls=':', color='k')                    # ligne verticale pointillee (':') noire ('k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')                                    
plt.title("L'ensemble de donnees")
plt.draw()


# resultat de classification 
plt.subplot(222)
plt.scatter( datax1, datax2, c = classe )
plt.quiver( 0, 0, 5*w[1]/np.linalg.norm(w), 5*w[2]/np.linalg.norm(w), angles='xy',scale_units='xy',scale=1 ) # vecteur normal
ax1, bx1 = -2, 8
# ax2, bx2 =  - w[0]*ax1 / w[1], - w[0]*bx1 / w[1]
ax2, bx2 =  - w[0] / w[2] - w[1]*ax1 / w[2] , - w[0] / w[2] - w[1]*bx1 / w[2]
plt.plot([ax1, bx1], [ ax2, bx2], 'g' )                 # hyperplan separateur
plt.axhline(0, ls=':', color='k')
plt.axvline(0, ls=':', color='k')
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')
plt.title('classification')
plt.draw()
      

print "# TRANSFORMATION D'ENTREES PAR LE RESEAU"
x1range  = np.arange(-10,10,0.1)
x2range  = np.arange(-10,10,0.1)
y = np.zeros(( len(x2range), len(x1range) ))
for i, x2 in enumerate(x2range):
    for j, x1 in enumerate(x1range):
        x       = [1, x1, x2]               # vecteur des entrees
        a       = np.dot( w, x)             # activation du neurone de sortie
        y[i,j]  = sigmoid(a)                # activite du neurone de sortie

# figure : transformation des entrees
X1, X2 = np.meshgrid(x1range, x2range)
ax = plt.subplot(223, projection='3d')
ax.plot_surface(X1, X2, y, rstride=10, cstride=10, cmap='jet')
plt.xlabel("x1")
plt.ylabel("x2")
ax.set_zlabel("y")
plt.title("Transformation d'entrees")
plt.draw()

plt.show()
exit()                                        
print "# FONCTION D'ERREUR"

w1range  = np.arange(-4,4,0.1)
w2range  = np.arange(-4,4,0.1)
err = np.zeros(( len(w2range), len(w1range), len(w0range) ))

for i, w2 in enumerate(w2range):
    for j, w1 in enumerate(w1range):
        w  = [w1, w2]
        e = 0
        for k in range(N):
            x = [ 1, datax1[k], datax2[k] ]
            a = np.dot( w, x)
            y = sigmoid(a)
            
            if labels[k] == 'r':
                target = 1
            else: 
                target = 0
                
        e = e + 0.5 * (y - target)**2
    
        err[i,j] = e


print "# figure : fonction d'erreur"
W1, W2 = np.meshgrid(w1range, w2range)
ax = plt.subplot(224, projection='3d')
ax.plot_surface(W1, W2, err, rstride=5, cstride=5, cmap='summer')
plt.xlabel("w1")
plt.ylabel("w2")
ax.set_zlabel("Erreur")
plt.title("Fonction-erreur")
plt.draw()

plt.show()
