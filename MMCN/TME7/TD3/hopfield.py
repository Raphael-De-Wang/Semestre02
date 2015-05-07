# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
from time import sleep

# nombre totale de neurones
N = 64

# fonction d'activation
#f = lambda r: 50*(1+np.tanh(r))
#f = lambda r: np.sign(r)
g = lambda a: 1./(1.+np.exp(-a))
    
# motif à mémoriser
#motif = np.zeros((8,8))
#motif[:,3:5] = 1
#motif[3:5, :] = 1

motif = np.array( [   [0, 0, 0, 1, 1, 0, 0, 0],
                      [0, 0, 0, 1, 1, 0, 0, 0],
                      [1, 1, 1, 1, 1, 1, 1, 1],
                      [0, 0, 0, 1, 1, 0, 0, 0],
                      [0, 0, 0, 1, 1, 0, 0, 0],
                      [0, 0, 1, 1, 1, 1, 0, 0],
                      [0, 1, 1, 1, 1, 1, 1, 0],
                      [1, 1, 1, 0, 0, 1, 1, 1]  ], dtype=float )
                      
# poids synpatiques hebbiens
eta = 0.01
W = eta*np.outer( motif.flatten(), motif.flatten() )

# pas d'entrées
x = 0

# PLOT
plt.figure(1); plt.clf(); plt.show()

# motif à mémoriser
plt.subplot(311)
plt.imshow(motif)
plt.draw()

# matrice de poids
plt.subplot(312); plt.cla()
plt.imshow(W)
plt.draw()

raw_input("\n Press any key ... \n")

# paramètres de simulation
dt = 0.1
time = np.arange(0,10,dt)
T = len(time)
y = np.zeros((T, N))        # pour chaque instant de tamps, il y a vecteur y de N valeurs

# condition initiale aléatoire
y[0] = np.random.rand(N)

for t in range(1, T):
    y[t] = y[t-1]+ dt*( -y[t-1] + g(np.dot(W, y[t-1]) + x) )

    print("Iteration ", t)
    
    plt.subplot(313); plt.cla()
    plt.imshow(y[t].reshape(8,8))
    plt.draw()
    
    sleep(0.1)

