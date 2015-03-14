# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt

# Temps

T   = 500.                  # ms, Temps total de la simulation
dt  = 1.                    # ms, Pas de temps de la simulation
t   = np.arange(0,T,dt)     # ms, Vecteur des temps en lesquels sont calculées les variables
N   = len(t)                # Longueur du vecteur des temps 

# Stimulation par un créneau de courant injecté

i_amp = 0.006   # uA, Courant injecté pendant le créneau
t_on  = 100     # ms, Début du créneau
t_off = 400     # ms, Fin du créneau

# Excitabilité

tau_m = 10      # ms, Constante de temps membranaire
R_m   = 4000    # KOhm, Résistance de membrane 
theta = 20      # mV, Seuil du potentiel d'action
v_P   = 100     # mV, Potentiel au pic du PA 

# Initialisation

i_inj  = np.zeros(N)
v      = np.zeros(N)

# Simulation

for i in range(1,N):

    # courant injecté
    if (t[i]>=t_on) and (t[i]<t_off):
        i_inj[i]  = i_amp
            
    # compléter la ligne suivante
    v[i] = v[i-1] + dt*(-v[i-1]+R_m*i_inj[i])/tau_m

    # detection de PA
    if v[i-1] > theta:
        v[i] = v_P

    # reset à 0
    if v[i-1] == v_P:
        v[i] = 0;


plt.figure(1); plt.show(); plt.clf();

# Figure Courant
plt.subplot(211)
plt.plot(t,i_inj);
plt.ylim(0, 0.010)
plt.xlabel('Temps (ms)')
plt.ylabel('Courant, (microA)')
plt.draw()


# Figure potentiel du neurone post-synaptique
plt.subplot(212)
plt.plot(t,v)
plt.ylim(0, 110)
plt.xlabel('Temps (ms)')
plt.ylabel('V (mV)')
plt.title('Potentiel membranaire post-synaptique')
plt.draw()
