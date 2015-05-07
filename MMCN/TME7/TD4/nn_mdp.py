# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

# le nombre d'epreuves
K = 1000

# le nombre de neurones d'entree (egal au nombre d'etats)
S = 10
D = 10

# le nombre de neurones de sortie (correspond au nombre d'actions)
N = 2

# parametres du modele
eta = 0.1
eps = 0.1
gam = 0.9

# matrice de poids
W = np.zeros((N,D), dtype=float)
# W = np.random.rand(N,D)
# W[:,-1] = [0,0]

Q = np.zeros((S,N))
mem_Q  = []


# boucle pour K epreuves
for i in range(K):

    print('Epreuve ' + str(i))

    # l'etat initial
    s = 0

    # l'etat est encode par l'activite du neurone d'entree correspondant
    x = np.zeros(D)
    x[s] = 1.

    # la fonction-valeur correspond a l'activite des neurones de sortie
    Q[s,:] = np.dot( W, x)

    # choix d'une action, strategie epsilon-greedy
    if np.random.rand() < eps:
        a = np.random.randint(N)    # action aleatoire, exploration
    else:
        a = np.argmax(Q[s,:])       # action optimale, exploitation

    # repeter jusqu'a la fin de l'epreuve
    while s!=S-1:

        # effectuer l'action choisie et passer a l'etat s`
        if a==0:    
            if s<S-1: s_new = s + 1     # aller a droite, mais pas depasser le mur
            else: s_new = s

        elif a==1:       
            if s>0: s_new = s - 1       # aller a gauche, mais pas depasser le mur
            else: s_new = s

        else:
            print("Action n'existe pas")

        # obtenir la recompense
        if s_new == 9: 
            r = 1.
        else:
            r = 0.

        # calculer l'activite des neurones d'entree et la valeur des actions dans le nouvel etat
        x_new = np.zeros(D)
        x_new[s_new] = 1.

        Q[s_new, :] = np.dot( W, x_new)

        # choisir la nouvelle action a` dans l'etat s`
        if np.random.rand() < eps :
            a_new = np.random.randint(N)    # action aleatoire, exploration
        else:
            a_new = np.argmax(Q[s_new,:])   # action optimale, exploitation

        # mettre a jour la matrice de poids
        delta = r + gam * Q[s_new, a_new] - Q[s,a]      # signal dopaminergique
        W[a, :] = W[a,:]  + eta*delta*x                 # plasticite synaptique

        # initialiser la nouvelle epreuve
        a = a_new
        s = s_new
        x = x_new

    # garder les resultats pour visualisation
    mem_Q.append(Q.copy())  


mem_Q = np.array(mem_Q)

# valeur de Q(s,a0) pendant l'apprentissage pour s=0, s=4 et s=8
Q_s0_a0 = mem_Q[:,0,0]
Q_s4_a0 = mem_Q[:,4,0]
Q_s8_a0 = mem_Q[:,8,0]

# valeur de Q(s,a0) apres l'apprentissage pour tous les etats et a=0
Q_s_a0 = mem_Q[-1,:,0]

# valeur theorique
Q_theor = gam**(8-np.arange(0,9))

# initialisation de graphisme
plt.figure(1); plt.clf()

plt.subplot(211)
plt.plot(Q_s0_a0, label = 's=0, a=0' )
plt.plot(Q_s4_a0, label = 's=4, a=0' )
plt.plot(Q_s8_a0, label = 's=8, a=0' )
plt.axhline(gam**8, ls = ':', lw=2, color='b')
plt.axhline(gam**4, ls = ':', lw=2, color='g')
plt.axhline(gam**0, ls = ':', lw=2, color='r')
plt.ylim(-0.1, 1.1)
plt.xlabel('epreuves')
plt.ylabel('Q(s,a)')
plt.legend()
plt.draw()

plt.subplot(212)
plt.plot(Q_s_a0[:-1], 'ko-', label = 'appris' )
plt.plot(Q_theor, 'ko:', label = 'theorique' )
plt.xlabel('etats')
plt.ylabel('Q(s,a0)')
plt.legend(loc='upper left')
plt.draw()

plt.show() 

