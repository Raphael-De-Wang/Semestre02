# -*- coding: utf-8 -*-

import numpy as np
from arftools import *


class Kmeans(object):
    def __init__(self,K):
        self.K = K
        self.RESET = True
        self.STOP  = False
    def fit(self, X, Y):
        self.data = X
        self.label= Y
        self.dim  = len(X[0])
        self.scale= np.max(X,axis=0)
        self.optimize()

    def init_prototype(self):
        self.prototype = np.random.randn(self.K,self.dim)*self.scale
        
    def optimize(self):
        # Initialisation des k Centres
        if self.RESET:
            self.init_prototype()

        self.iCount = 0
        while not self.STOP:
            self.iCount = self.iCount + 1
            # Affectation des points
            C = self.predict(self.data)
            # Estimation des centres
            self.update_prototype(C)
            print self.prototype
            if self.iCount > 10:
                self.STOP = True

    def update_prototype(self,affection):
        try:
            self.prototype = np.array([ self.class_cost(i,affection==i) for i,mu in enumerate(self.prototype)])
        except:
            print i, affection
        
    def score(self):
        pass

    def predict(self,X):
        return np.argmin([ [ np.linalg.norm(x-mu)**2 for mu in self.prototype ] for x in X ], axis=1)

    def class_cost(self,class_indice, data_indice):
        mu = self.prototype[class_indice]
        return np.sum([ [ np.linalg.norm(x[i]-mu[i])**2 for i in range(self.dim) ] for x in self.data[data_indice]], axis=0)/sum(data_indice)
    
    def cost(self):
        return np.sum([ [ np.linalg.norm(x-mu) for mu in self.prototype ] for x in self.data ])
        
t = 3
X,Y = gen_arti(data_type=t)

cluster = Kmeans(2)
cluster.fit(X,Y)
