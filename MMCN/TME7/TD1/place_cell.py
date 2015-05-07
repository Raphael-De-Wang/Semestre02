# -*-coding:Latin-1 -*
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# environnement
LX = 1000        # mm
LY = 1000        # mm

# cellule de lieu
A = 100.         # Hz, activite maximale
S = 100.         # mm, largeur du champs recepteur
x_pref = 400.    # mm, coordonnee x de la position preferee
y_pref = 600.    # mm, coordonnee y de la position preferee

posx = np.arange(0, LX, 20)
posy = np.arange(0, LY, 20)
NX = len(posx)
NY = len(posy)

# enregistrement d'une cellule de lieu dans l'hippocampe
r = np.zeros((NY,NX))

for i in range(1,NY):
    for j in range(1,NX):
        
        # position du rat dans l'environnement
        x = posx[j]
        y = posy[i]

        # activite d'une cellule de lieu dans la position (x,y)
        r[i,j] = A*np.exp( -(x_pref-x)**2/(2*S**2) - (y-y_pref)**2/(2*S**2) )

# Figure 
plt.figure(1); plt.show(); plt.clf()
 
plt.subplot(211)
plt.contourf(posy, posx, r)
plt.axis('scaled')
plt.xlim(0,LX)
plt.ylim(0,LY)
plt.title('Activite d`UNE CELLULE de lieu')
plt.draw()

ax = plt.subplot(212, projection = '3d')
[X,Y] = np.meshgrid(posx, posy)
ax.plot_surface(X, Y, r, rstride=1, cstride=1, cmap='jet')
plt.title('Activite d`UNE CELLULE de lieu')
plt.draw()

