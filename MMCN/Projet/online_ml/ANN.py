# -*- coding: utf-8 -*-
#!env python

import numpy as np
import matplotlib.pyplot as plt
from mytools import *

np.set_printoptions(precision=4)
# np.seterr(all='raise')


def sigmoid(h):
    try:
        return 1./(1.+np.exp(-h))
    except:
        print h
        raise

def grad_sigmoid(h):
    return sigmoid(h)*(1.-sigmoid(h))

class Neuron(object):
    def __init__(self,inputSize,f=None,grad_f=None):
        # print "new neuron."
        self._f        = f
        self._grad_f   = grad_f
        self.inputSize = inputSize
        self.weights   = self.x_random()
    def x_random(self):
        return np.random.random(self.inputSize)-0.5
    def a(self,x):
        self.A = np.sum(self.weights*x)
        return self.A
    def f(self,x):
        try:
            return self._f(self.a(x))
        except:
            print self.weights
            raise
    def grad_f(self,x):
        return self._grad_f(self.a(x))
    def updateWeights(self,delta,prevZ,eta=0.01):
        try:
            self.weights = self.weights - delta * prevZ * eta
        except:
            print detla
            raise
    
class Layer(object):
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.01):
        self.g         = g
        self.grad_g    = grad_g
        self.inputSize = inputSize
        self.layerSize = layerSize
        self.eta       = eta
        self.delta     = np.array([np.inf])
    def forward(self,prevZ):
        self.prevZ  = prevZ
        self.Z      = np.array([ n.f(prevZ) for n in self.layer ])
        self.grad_Z = np.array([ n.grad_f(prevZ) for n in self.layer ])
        return self.Z
    def getWeightMatrix(self):
        return np.array([ n.weights for n in self.layer ])
    def updateWeights(self):
        for i in range(self.layerSize):
            self.layer[i].updateWeights(self.delta[i],self.prevZ,self.eta)
    def calculateDelta(self,nextDeltaWeightSum):
        raise NotImplementedError("calculateDelta non implemente")


class HiddenLayer(Layer):
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.01):
        super(HiddenLayer,self).__init__(inputSize,layerSize,g,grad_g)
        self.layer = np.array([ Neuron(inputSize,g,grad_g) for i in xrange(layerSize-1)])
        bias       = Neuron(inputSize,lambda x: 1 ,lambda x: 0 )
        self.layer = np.append(self.layer,bias)
    def calculateDelta(self,nextDeltaWeightSum):
        self.delta = self.grad_Z * nextDeltaWeightSum
        # print self.delta
        return self.delta

class OutputLayer(Layer):
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.01):
        super(OutputLayer,self).__init__(inputSize,layerSize,g,grad_g)
        self.layer     = np.array([ Neuron(inputSize,g,grad_g) for i in xrange(layerSize)])
    def calculateDelta(self,y):
        self.delta = np.array([yk - ykCible(self.x,self.r,dk(i,len(self.Z))) for i,yk in enumerate(self.Z)])
        # print self.delta
        return self.delta
        
class Classifier(object):
    def fit(self,x,y):
        raise NotImplementedError("fit non  implemente")
    def predict(self,x):
        raise NotImplementedError("predict non implemente")
    def score(self,x,y):
        return (self.predict(x)==y).mean()

class NeuroNetwork(Classifier):
    def __init__(self,topology,K,eta=0.01):
        self.err = []
        self.topology    = topology
        self.hiddenLayers= np.array([ HiddenLayer(inputSize,layerSize,g,grad_g,eta) for g,grad_g,inputSize,layerSize in self.topology])
        self.outputLayer = OutputLayer(layerSize,K,lambda x: x,lambda x: 1,eta)
        self.layers      = np.append(self.hiddenLayers,self.outputLayer)
        self.eta         = eta
    def fit(self,x,y):
        self.x   = x 
        self.y   = y
        self.optimize()
    def forward(self,X):
        Z = input_layer_active(X)
        Z = np.append(Z,1)
        for l in self.layers:
            Z = l.forward(Z)
        self.layers[-1].x = X[0]
        self.layers[-1].r = X[-1]
        return Z
    def backProp(self,x,y):
        deltaWeiSum = y
        # Layers
        for i in range(len(self.layers)-1,-1,-1):
            delta       = self.layers[i].calculateDelta(deltaWeiSum)
            wm          = self.layers[i].getWeightMatrix()
            # deltaWeiSum = np.sum(delta * wm.T,axis=1)
            deltaWeiSum = np.dot(wm.T,delta)
        # update layer weights
        for i in range(len(self.layers)):
            self.layers[i].updateWeights()
    def optimize(self):
        for i in range(len(self.x)):
            self.forward(self.x[i])
            self.backProp(self.x[i],self.y[i])
            self.err.append(sum(self.outputLayer.delta**2))
            if self.err[-1] < 1e-2:
                return
    def predict(self,X):
        yk = []
        ykCible = []
        for x in X:
            self.forward(x)
            # print self.outputLayer.Z
            yk.append(2.*np.argmax(self.outputLayer.Z)-90)
            ykCible.append(x[0]+x[-1])
        return np.array(yk),np.array(ykCible)
    def error(self):
        return abs(self.outputLayer.delta.sum())
