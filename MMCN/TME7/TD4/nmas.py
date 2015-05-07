# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

# le nombre d'epreuves
K = 5000

# le nombre d'etats 
# le nombre de machines a sous, chacune correspond a une couleur du crois de fixation
S = 3

# le nombre d'actions (direction d'une saccade) dans chaque etat 
N = 4

# recompenses moyennes (le nombre de gouttes de jus)
E_r = np.random.randint(0,20, (S,N))

# parametres du modele
eta = 0.1
eps = 0.1

# initialisation
Q = np.random.rand(S,N)   # fonction valeur

mem_Q  = []

# boucle pour K epreuves
for i in range(K):

    # choix aleatoire d'une machine a sous (ou d'un etat)
    s = np.random.randint(S)

    # choix d'une action, strategie epsilon-greedy
    if np.random.rand() < eps :
        a = np.random.randint(N)    # action aleatoire, exploration
    else:
        a = np.argmax(Q[s,:])       # action optimale, exploitation

    # recompense
    r =  np.random.poisson(E_r[s,a])

    # mettre a jour l'estimation de la fonction valeur
    delta = r - Q[s,a]
    Q[s,a] = Q[s,a] + eta*delta

    # garder les resultats pour visualisation
    mem_Q.append(Q.copy())  


# initialisation de graphisme
plt.figure(1); plt.clf()

# estimation de la fonction valeur

mem_Q = np.array(mem_Q)

Q_mas0 = mem_Q[:, 0, :]
Q_mas1 = mem_Q[:, 1, :]
Q_mas2 = mem_Q[:, 2, :]

plt.subplot(311)
plt.plot(Q_mas0 )
plt.axhline(E_r[0,0], ls = ':', lw=2, color='b')
plt.axhline(E_r[0,1], ls = ':', lw=2, color='g')
plt.axhline(E_r[0,2], ls = ':', lw=2, color='r')
plt.axhline(E_r[0,3], ls = ':', lw=2, color='c')
plt.xlabel('epreuves')
plt.ylabel('Q(a)')
plt.title('Machine a sous 1')
plt.draw()

plt.subplot(312)
plt.plot(Q_mas1 )
plt.axhline(E_r[1,0], ls = ':', lw=2, color='b')
plt.axhline(E_r[1,1], ls = ':', lw=2, color='g')
plt.axhline(E_r[1,2], ls = ':', lw=2, color='r')
plt.axhline(E_r[1,3], ls = ':', lw=2, color='c')
plt.xlabel('epreuves')
plt.ylabel('Q(a)')
plt.title('Machine a sous 2')
plt.draw()


plt.subplot(313)
plt.plot(Q_mas2 )
plt.axhline(E_r[2,0], ls = ':', lw=2, color='b')
plt.axhline(E_r[2,1], ls = ':', lw=2, color='g')
plt.axhline(E_r[2,2], ls = ':', lw=2, color='r')
plt.axhline(E_r[2,3], ls = ':', lw=2, color='c')
plt.xlabel('epreuves')
plt.ylabel('Q(a)')
plt.title('Machine a sous 3')
plt.draw()

plt.show()
