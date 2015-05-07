# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

# poids synaptiques récurrents et les entrées
w = -10
x = 5.0

# fonction d'activation
g = lambda a: 1./(1+np.exp(-a))

# écrire les nullclines du système
ncy1 = lambda y: -(x + np.log((1-y)/y))/w
ncy2 = lambda y: g(w*y + x)

# paramètres de la simulation
dt = 0.1
time = np.arange(0,10,dt)
T = len(time)
y = np.c_[np.zeros(T), np.zeros(T)]

# Quelle est la matrice des poids synaptiques ?
W = np.array([[0,w],[w,0]])

# condition initiale
y[0] = [0.2, 0.4]

for t in range(1, T):

    # méthode d'Euler pour dy/dt=-y+g(wy+x)
    y[t] = y[t-1]+ dt*(-y[t-1] + g( np.dot(W, y[t-1])+x))


# initialisation de graphisme
plt.figure(1); plt.show(); plt.clf();

y1_interv1 = np.arange( -0.1, 1.3, 0.01)
y1_interv2 = np.arange( 0.001, 0.999, 0.001)

plt.subplot(211)
plt.plot(y1_interv2, ncy1(y1_interv2), 'g-', lw=3, label="nullcline y1'" )
plt.plot(y1_interv1, ncy2(y1_interv1), 'b-', lw=3, label="nullcline y2'" )
plt.plot( y[:,0], y[:,1], 'k-', lw=3)
plt.plot( y[0,0], y[0,1], 'ro', lw=3)
plt.axis('scaled')
plt.axhline(0, ls='-', color = 'k')
plt.axvline(0, ls='-', color = 'k')
plt.legend()
plt.draw()


plt.subplot(212); 
plt.plot(time, y[:,0], 'k-' , label='y1(t)')
plt.plot(time, y[:,1], 'k--', label='y2(t)' )
plt.axhline(0, ls=':', color = 'k')
plt.axvline(0, ls=':', color = 'k')
plt.legend()
plt.draw()
