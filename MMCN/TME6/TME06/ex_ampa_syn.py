# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

T   = 30.                   # ms, Temps total de la simulation
dt  = 0.1                   # ms, Pas de temps de la simulation
t   = np.arange(0.,T,dt)    # ms, Vecteur des temps en lesquels sont calculées les variables
N   = len(t)                # Longueur du vecteur des temps 

# Stimulation pré-synaptique

t_pre        = 5.                    # ms, Temps d'arrivée du PA pré-synaptique
ind_pre      = np.round(t_pre/dt)    # Indice correspondant
pre          = np.zeros(N)           # Vecteur des PA pré-synaptiques (0=rien 1=PA)
pre[ind_pre] = 1                     # Un PA au temps t_pre (indice ind_pre)

# diffusion de neurotransmetteur

T_amp = 1.          # mM, concentration maximale du neurotransmetteur
T_dur = 1.          # ms, 1ms / durée de l'impulsion de du neurotransmetteur

# Récepteurs AMPA

alpha_ampa = 0.93     # 1/(mM x ms)
beta_ampa  = 0.19     # 1/ms

# Initialisation

T_nt = np.zeros(N)
p_ampa = np.zeros(N)

# Simulation

for i in range(1,N):
    
    # neurotransmetteur
    if (t[i]>=t_pre) and (t[i]<t_pre+T_dur):
        T_nt[i] = T_amp
    
    # ajouter le modèle des récepteurs AMPA
    
    
plt.figure(1); plt.show(); plt.clf()

# Figure PA pre-synaptique
plt.subplot(311)
plt.stem(t, pre, markerfmt='.')
plt.ylim(-0.1, 5)
plt.xlabel('Temps (ms)')
plt.ylabel('PA pre-synaptique')
plt.draw()

# Figure neurotransmetteur
plt.subplot(312)
plt.plot(t,T_nt)
plt.ylim(0, 2)
plt.xlabel('Temps (ms)')
plt.ylabel('concentration du neurotransmetteur')
plt.draw()

# Figure p(ouverture) du canal AMPA
plt.subplot(313)
plt.plot(t,p_ampa)
plt.ylim(0, 2)
plt.xlabel('Temps (ms)')
plt.ylabel('p(ouverture) AMPA')
plt.draw()
