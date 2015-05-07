# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

# données
# donnees = np.array([[ 1., 1.],  [ 3., 1.], [-2., 0.], [-2.,-2.]]) 
# N = len(donnees)

# générer des données avec una matrice de covariance donnée
N = 50
moyenne = np.array([0.0, 0.0])
C = np.array([ [ 1, -0.7], [ -0.7, 1] ])
donnees = np.random.multivariate_normal(moyenne, C, size=N)

cov = np.cov(donnees.T)

x1 = donnees[:,0]
x2 = donnees[:,1]
'''
# calculer la matrice de covariance
cov_x1x1 = cov[0,0]
cov_x1x2 = cov[1,0]
cov_x2x2 = cov[1,1]
'''

# C = np.array([[cov_x1x1, cov_x1x2], [cov_x1x2, cov_x2x2]])
C = cov

# valeurs et vecteurs propres
valp, vecp =  np.linalg.eig(C)
print('Valeurs propres de la matrice de covariance: %.1f, %.1f' % (valp[0], valp[1]) )

# l'axe principal
ind_max = np.argmax( valp)      # l'index de la valeur propre maximale
axe_princ = vecp[:, ind_max]    # le vecteur propre correspondant


plt.figure(1); plt.clf(); plt.show()

plt.subplot(211)
plt.scatter( x1, x2 )  
plt.quiver( 0, 0, axe_princ[0], axe_princ[1], angles='xy',scale_units='xy',scale=1, color='r' ) # axe principal
plt.axhline(0, ls=':', color='k')                    
plt.axvline(0, ls=':', color='k')                    
plt.xlabel("x1")
plt.ylabel("x2")
plt.axis('scaled')                                    
plt.draw()

raw_input("\n Press any key \n")

# réseau avec l'apprentissage Hebbien
T = 100
eta = 0.05
w = np.random.rand(2) - 0.5

for i in range(T):
    
    p = np.random.randint(N)
    x = donnees[p]
    y = np.dot(w, x)

    # règle d'apprentissage hebbien
    # w = w + eta*y*x
    w = w + eta*y*(x-y*w)

    print('Iteration ', i)
    plt.subplot(212); #plt.cla()
    plt.scatter( x1, x2 )  
    plt.plot([0, w[0]], [0, w[1]], '-r')
    plt.axhline(0, ls=':', color='k')                    
    plt.axvline(0, ls=':', color='k')                    
    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.axis('scaled')                                    
    plt.draw()

    sleep(0.1)

plt.subplot(211)
plt.quiver( 0, 0, w[0], w[1], angles='xy',scale_units='xy',scale=1, color='g' ) # axe principal
plt.draw()

