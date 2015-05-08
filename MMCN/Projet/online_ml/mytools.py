# -*- coding: utf-8 -*-
import numpy as np

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

toInt = lambda x : np.floor(x) if x < 0 else np.ceil(x)

a_unif= lambda i,N1 : 45*(i-(N1-1)/2.)/((N1-1)/2.)
b_unif= lambda i,N2 : i/(N2-1.)
c_unif= lambda i,N3 : i/(N3-1.)

sigma = 10
posAct= lambda x,indice,N1 : np.exp(-(x-a_unif(indice,N1)**2)/(2.*(sigma**2)))
rPlAct= lambda r,indice,N2 : b_unif(indice,N2) * ( r + 45. ) / 90.
rMiAct= lambda r,indice,N3 : c_unif(indice,N3) * ( 45. - r ) / 90.

ykCible = lambda x,r,dk : np.exp(-((x+r-dk)**2)/(2.*sigma**2))

def data_uniform(N1=41,N2=11,N3=11,OUTPUT=91):
    position = np.array([ p for p in range(-45,46) for repeat in range(-45,46) ])*1.
    r        = np.array([ r for repeat in range(-45,46) for r in range(-45,46) ])*1.
    label    = np.zeros((len(r),OUTPUT))
    for i,l in enumerate(map(toInt,(position + r)/2)):
        label[i,l] = 1 + N1/2
    pos    = np.ones((N1,len(r))) * position
    r_plus = np.ones((N2,len(r))) * r
    r_minus= np.ones((N3,len(r))) * r
    data   = np.vstack((pos,r_plus,r_minus))
    return data.T, label

def data_random(N1=41,N2=11,N3=11,OUTPUT=91,data_size=10000):
    position = np.random.randint(-45,45,data_size)*1.
    r        = np.random.randint(-45,45,data_size)*1.
    # label  = map(toInt,(position + r)/2)
    label    = np.zeros((len(r),OUTPUT))
    for i,l in enumerate(map(toInt,(position + r)/2)):
        label[i,l] = 1 + N1/2
    pos    = np.ones((N1,data_size)) * position
    r_plus = np.ones((N2,data_size)) * r
    r_minus= np.ones((N3,data_size)) * r
    data   = np.vstack((pos,r_plus,r_minus))
    return data.T, label

def input_layer_active(D,N1=41,N2=11,N3=11):
    d = D.T
    ones = np.ones(len(d[0]))
    for i,data in enumerate(d):
        if i < N1:
            d[i] = map(posAct,data,ones*i,ones*N1)
        elif i < N1 + N2:
            d[i] = map(rPlAct,data,ones*i,ones*N2)
        else:
            d[i] = map(rMiAct,data,ones*i,ones*N3)
    return d.T
