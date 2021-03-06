# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

# le nombre d'epreuves
K = 2000

# le nombre d'actions (direction d'une saccade)
N = 4

# recompenses moyennes (le nombre de gouttes de jus)
E_r = np.array([ 2., 4., 6., 8.])

# parametres du modele
eta = 0.1
eps = 0.1

# initialisation
Q = np.random.rand(N)   # fonction valeur

mem_Q  = []
mem_delta = []

# boucle pour K epreuves
for i in range(K):

    # choix d'une action, strategie epsilon-greedy
    if np.random.rand() < eps :
        a = np.random.randint(N)    # action aleatoire, exploration
    else:
        a = np.argmax(Q)            # action optimale, exploitation

    # recompense
    r =  np.random.poisson(E_r[a])

    # mettre a jour l'estimation de la fonction valeur
    delta = r - Q[a]
    Q[a] = Q[a] + eta*delta

    # garder les resultats pour visualisation
    mem_Q.append(Q.copy())  

    if a == 2:   # saccade vers la premiere cible
        mem_delta.append(delta)


# initialisation de graphisme
plt.figure(1); plt.clf()

# distribution de recompense
r0 = np.random.poisson(E_r[0], 1000)
r1 = np.random.poisson(E_r[1], 1000)
r2 = np.random.poisson(E_r[2], 1000)
r3 = np.random.poisson(E_r[3], 1000)

plt.subplot(321)
plt.hist(r0, color = 'b', label = 'r1')
plt.axvline(E_r[0], lw=2, color='k')
plt.legend()
plt.draw()

plt.subplot(322)
plt.hist(r1, color = 'g', label = 'r2')
plt.axvline(E_r[1], lw=2, color='k')
plt.legend()
plt.draw()

plt.subplot(323)
plt.hist(r2, color = 'r', label = 'r3')
plt.axvline(E_r[2], lw=2, color='k')
plt.legend()
plt.draw()

plt.subplot(324)
plt.hist(r3, color = 'c', label = 'r4')
plt.axvline(E_r[3], lw=2, color='k')
plt.legend()
plt.draw()


# estimation de la fonction valeur
plt.subplot(325)
plt.plot( np.array(mem_Q))
plt.axhline(E_r[0], ls = ':', lw=2, color='b')
plt.axhline(E_r[1], ls = ':', lw=2, color='g')
plt.axhline(E_r[2], ls = ':', lw=2, color='r')
plt.axhline(E_r[3], ls = ':', lw=2, color='c')
plt.xlabel('epreuves')
plt.ylabel('Q(a)')
plt.draw()

# erreur de prevision
plt.subplot(326)
plt.stem( range(len(mem_delta)),mem_delta)
plt.xlabel('epreuves quand a=2 est choisie')
plt.axhline(0, ls=':', color='k')
plt.ylabel('delta')
plt.draw()

plt.show()
