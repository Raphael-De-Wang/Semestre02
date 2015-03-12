# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
    
T   = 1000.                 # ms, Temps total de la simulation
dt  = 1.                    # ms, Pas de temps de la simulation
t   = np.arange(0,T,dt)     # ms, Vecteur des temps en lesquels sont calculées les variables
N   = len(t)             # Longueur du vecteur des temps 

# Stimulation pré-synaptique

f_pre       = 100.             # Hz, Fréquence des PA pré-synaptiques
p_pre       = f_pre*dt/1000    # Probabilité d'avoir un PA dans un pas de temps
pre         = np.random.rand(N)<p_pre  # Vecteur des PA pré-synaptiques Poisson (0=rien 1=PA)

# Canaux récepteurs AMPA

tau_ampa = 5.26   # ms, Constante de temps de fermeture des canaux AMPA
p_max    = 1      # Incrément d'ouverture des canaux lors d'un PA pré-synaptique
g_ampa   = 0.1    # mS/cm2, conductance maximale de recepteurs AMPA
E_ampa   = 60     # mV, conductance maximale de recepteurs AMPA

# Excitabilité

tau_m =  10.  # ms, Constante de temps membranaire
r_m   =  10.  # KOhm.cm2, Résistance de fuite 
theta =  20.  # mV, Seuil du potentiel d'action
v_P   =  100.  # mV, Potentiel au pic du PA 

# Initialisation

i_ampa = np.zeros(N)
p_ampa = np.zeros(N)
v      = np.zeros(N)

# Simulation

for i in range(1,N):

    p_ampa[i] = p_ampa[i-1] - dt / tau_ampa * p_ampa[i-1] 
    
    if pre[i] == 1:
        p_ampa[i] = p_max

    # courant AMPA
    i_ampa[i] = g_ampa * p_ampa[i] * (v[i-1] - E_ampa) 

    # potentiel membranaire
    v[i] = v[i-1] + dt / tau_m * ( - v[i-1] - r_m * i_ampa[i] )

    # detection de PA
    if v[i-1] > theta:
        v[i] = v_P

    # reset à 0
    if v[i-1] == v_P:
        v[i] = 0

plt.figure(1); plt.show(); plt.clf()

# Figure PA pre-synaptique
plt.subplot(311)
plt.stem(t, pre, markerfmt='.' )
plt.ylim( -0.1, 2)
plt.xlabel('Temps (ms)')
plt.ylabel('PA pre-synaptique')
plt.draw()


# Figure p(ouverture) du canal AMPA
plt.subplot(312)
plt.plot(t,p_ampa)
plt.ylim(-0.1, 2)
plt.xlabel('Temps (ms)')
plt.ylabel('p(ouverture) AMPA')
plt.draw()

# Figure potentiel du neurone post-synaptique
plt.subplot(313)
plt.plot(t,v)
plt.axhline(theta, ls=':', color='k')
plt.ylim(0, 110)
plt.xlabel('Temps (ms)')
plt.ylabel('V (mV)')
plt.draw()
