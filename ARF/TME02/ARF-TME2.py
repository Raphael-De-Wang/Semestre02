# -*- coding: utf-8 -*-
import matplotlib
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from collections import Counter

def gen_arti(centerx=1,centery=1,sigma=0.1,nbex=1000,data_type=0,epsilon=0.02):
    #center : entre des gaussiennes
    #sigma : ecart type des gaussiennes
    #nbex : nombre d'exemples
    # ex_type : vrai pour gaussiennes, faux pour echiquier
    #epsilon : bruit

    if data_type==0:
        #melange de 2 gaussiennes
        xpos=np.random.multivariate_normal([centerx,centery],np.diag([sigma,sigma]),nbex/2)
        xneg=np.random.multivariate_normal([-centerx,-centery],np.diag([sigma,sigma]),nbex/2)
        data=np.vstack((xpos,xneg))
        y=np.hstack((np.ones(nbex/2),-np.ones(nbex/2)))
        
    if data_type==1:
        #melange de 4 gaussiennes
        xpos=np.vstack((np.random.multivariate_normal([centerx,centery],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([-centerx,-centery],np.diag([sigma,sigma]),nbex/4)))
        xneg=np.vstack((np.random.multivariate_normal([-centerx,centery],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([centerx,-centery],np.diag([sigma,sigma]),nbex/4)))
        data=np.vstack((xpos,xneg))
        y=np.hstack((np.ones(nbex/2),-np.ones(nbex/2)))

    if data_type==2:
        #echiquier
        data=np.reshape(np.random.uniform(-4,4,2*nbex),(nbex,2))
        y=np.ceil(data[:,0])+np.ceil(data[:,1])
        y=2*(y % 2)-1

    # un peu de bruit
    data[:,0]+=np.random.normal(0,epsilon,nbex)
    data[:,1]+=np.random.normal(0,epsilon,nbex)
    # on mélange les données
    idx = np.random.permutation((range(y.size)))
    data=data[idx,:]
    y=y[idx]
    return data,y

data,y = gen_arti()
    
#affichage en 2D des donnees
def plot_data(x,labels):
    fig = plt.figure()
    plt.scatter(x[labels<0,0],x[labels<0,1],c='red',marker='x')
    plt.scatter(x[labels>0,0],x[labels>0,1],c='green',marker='+')
    plt.show()
    
plot_data(data,y)

#Frontiere de decision
def plot_frontiere(x,f,step=20): # script qui engendre une grille sur l'espace des exemples, calcule pour chaque point le label
                                # et trace la frontiere
    mmax=x.max(0)
    mmin=x.min(0)
    x1grid,x2grid=np.meshgrid(np.linspace(mmin[0],mmax[0],step),np.linspace(mmin[1],mmax[1],step))
    grid=np.hstack((x1grid.reshape(x1grid.size,1),x2grid.reshape(x2grid.size,1)))
    # calcul de la prediction pour chaque point de la grille
    res=np.array([f(grid[i,:]) for i in range(x1grid.size)])
    res=res.reshape(x1grid.shape)
    # tracer des frontieres
    plt.contourf(x1grid,x2grid,res,colors=('gray','blue'),levels=[-1,0,1])

class Classifier(object):
    """ Classe generique d'un classifieur
        Dispose de 3 méthodes :
            fit pour apprendre
            predict pour predire
            score pour evaluer la precision
    """

    def fit(self,x,y):
        raise NotImplementedError("fit non implemente")
    def predict(self,x):
        raise NotImplementedError("predict non implemente")
    def score(self,x,y):
        return (self.predict(x) == y).mean()
    
class Knn(Classifier):
    def __init__(self, k):
        self.k = k
        
    def fit(self,trX,trY):
        self.trX = trX
        self.trY = trY
        
    def predict(self,X):
        
            

class Parzen(Classifier):
    def __init__(self, noyau):
        self.noyau = noyau
        
    def fit(self,X,Y):
        self.X = X
        self.Y = Y
        
    def predict(self,X):
        for x in X:
            pass
