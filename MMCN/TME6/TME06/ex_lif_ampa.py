# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

T   = 50.                   # ms, Temps total de la simulation
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

alpha_ampa = 0.93   # 1/(mM x ms)
beta_ampa  = 0.19   # 1/ms
g_ampa     = 0.1    # mS/cm2, conductance maximale de recepteurs AMPA
# E_ampa     = 60.    # mV, Potentiel d'inversion du récépteur AMPA
E_ampa     = -10.    # mV, Potentiel d'inversion du récépteur AMPA

# Excitabilité

tau_m =  10.        # ms, Constante de temps membranaire
r_m   =  5.         # KOhm.cm2, Résistance de fuite 
theta =  20.        # mV, Seuil du potentiel d'action
v_P   =  100.       # mV, Potentiel au pic du PA 

# Initialisation

T_nt   = np.zeros(N)
p_ampa = np.zeros(N)
i_ampa = np.zeros(N)
v      = np.zeros(N)

# Simulation

for i in range(1,N):

    # neurotransmetteur
    if (t[i]>=t_pre) and (t[i]<t_pre+T_dur):
        T_nt[i] = T_amp
    
    # ajouter les recepterus AMPA
    tau_p = 1./(alpha_ampa*T_nt[i]+beta_ampa)
    p_inf = alpha_ampa*T_nt[i]/(alpha_ampa*T_nt[i]+beta_ampa)
    p_ampa[i] = p_ampa[i-1] + (p_inf-p_ampa[i-1])*dt/tau_p

    # compléter la ligne avec le calcul de courant AMPA
    i_ampa[i] = g_ampa*p_ampa[i]*(v[i]-E_ampa)

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
plt.xlabel('Temps (ms)')
plt.ylabel('V (mV)')
plt.draw()
