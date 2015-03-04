# -*- coding: utf-8 -*-

import numpy as np
from numpy import random
import matplotlib.pyplot as plt



def to_array(x):
    """ convert an vector to array if needed """
    if len(x.shape)==1:
        x=x.reshape(1,x.shape[0])
    return x

def gen_arti(centerx=1,centery=1,sigma=0.1,nbex=1000,data_type=0,epsilon=0.02):
         #center : entre des gaussiennes
         #sigma : ecart type des gaussiennes
         #nbex : nombre d'exemples
         # ex_type : vrai pour gaussiennes, faux pour echiquier
         #epsilon : bruit

         if data_type==0:
             #melange de 2 gaussiennes
             xpos=np.random.multivariate_normal([centerx,centerx],np.diag([sigma,sigma]),nbex/2)
             xneg=np.random.multivariate_normal([-centerx,-centerx],np.diag([sigma,sigma]),nbex/2)
             data=np.vstack((xpos,xneg))
             y=np.hstack((np.ones(nbex/2),-np.ones(nbex/2)))
         if data_type==1:
             #melange de 4 gaussiennes
             xpos=np.vstack((np.random.multivariate_normal([centerx,centerx],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([-centerx,-centerx],np.diag([sigma,sigma]),nbex/4)))
             xneg=np.vstack((np.random.multivariate_normal([-centerx,centerx],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([centerx,-centerx],np.diag([sigma,sigma]),nbex/4)))
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

    #affichage en 2D des donnees
def plot_data(x,labels):
        plt.scatter(x[labels<0,0],x[labels<0,1],c='red',marker='x')
        plt.scatter(x[labels>0,0],x[labels>0,1],c='green',marker='+')

def make_grid(xmin=-5,xmax=5,ymin=-5,ymax=5,data=None,step=20):
    if data!=None:
        xmax=np.max(data[:,0])
        xmin=np.min(data[:,0])
        ymax=np.max(data[:,1])
        ymin=np.min(data[:,1])
    x=np.arange(xmin,xmax,(xmax-xmin)*1./step)
    y=np.arange(ymin,ymax,(ymax-ymin)*1./step)
    xx,yy=np.meshgrid(x,y)
    grid=np.c_[xx.ravel(),yy.ravel()]
    return grid,xx,yy

    #Frontiere de decision
def plot_frontiere(x,f,step=20): # script qui engendre une grille sur l'espace des exemples, calcule pour chaque point le label
                                    # et trace la frontiere
        grid,xvec,yvec=make_grid(data=x,step=step)
        #mmax=x.max(0)
        #mmin=x.min(0)
        #x1grid,x2grid=np.meshgrid(np.linspace(mmin[0],mmax[0],step),np.linspace(mmin[1],mmax[1],step))
        #grid=np.hstack((x1grid.reshape(x1grid.size,1),x2grid.reshape(x2grid.size,1)))
        # calcul de la prediction pour chaque point de la grille
        res=f(grid)
        res=res.reshape(xvec.shape)
        # tracer des frontieres
        plt.contourf(xvec,yvec,res,colors=('gray','blue'),levels=[-1,0,1])
        
##################################################################

class Classifier(object):
    """ Classe generique d'un classifieur
        Dispose de 3 méthodes :
            fit pour apprendre
            predict pour predire
            score pour evaluer la precision
    """

    def fit(self,x,y):
        raise NotImplementedError("fit non  implemente")
    def predict(self,x):
        raise NotImplementedError("predict non implemente")
    def score(self,x,y):
        return (self.predict(x)==y).mean()
