# -*- coding: utf-8 -*-
#!env python

import numpy as np
import matplotlib.pyplot as plt

from arftools import *

def sigmoid(h):
    return 1./(1.+np.exp(-h))
def grad_sigmoid(h):
    return sigmoid(h)*(1.-sigmoid(h))

def tangH(h):
    return np.tanh(h)
def grad_tangH(h):
    return 1. - tangH(h)**2
    
class Neuron(object):
    def __init__(self,inputSize,f=None,grad_f=None):
        print "new neuron."
        self._f        = f
        self._grad_f   = grad_f
        self.inputSize = inputSize
        self.weights   = self.x_random()
    def x_random(self,low=-1,high=1):
        return np.random.randn(self.inputSize)*(high-low)+low
    def a(self,x):
        self.A = np.sum(self.weights*x)
        return self.A
    def f(self,x):
        return self._f(self.a(x))
    def grad_f(self,x):
        return self._grad_f(self.a(x))
    def updateWeights(self,delta,prevZ,eta=0.1,alpha=1):
        # self.weights = self.weights - self.x * delta * eta
        self.weights = self.weights - delta * prevZ * eta
        # w = w - eps*dcdw
    
class Layer(object):
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.2):
        self.g         = g
        self.grad_g    = grad_g
        self.inputSize = inputSize
        self.layerSize = layerSize
        self.eta       = eta
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
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.2):
        super(HiddenLayer,self).__init__(inputSize,layerSize,g,grad_g)
        self.layer     = np.array([ Neuron(inputSize,g,grad_g) for i in xrange(layerSize-1)])
        bias           = Neuron(inputSize,lambda x: 1 ,lambda x: 0 )
        self.layer     = np.append(self.layer,bias)
    def calculateDelta(self,nextDeltaWeightSum):
        self.delta = self.grad_Z * nextDeltaWeightSum
        return self.delta
        
class OutputLayer(Layer):
    def __init__(self,inputSize,layerSize,g,grad_g,eta=0.2):
        super(OutputLayer,self).__init__(inputSize,layerSize,g,grad_g)
        self.layer     = np.array([ Neuron(inputSize,g,grad_g) for i in xrange(layerSize)])
    def calculateDelta(self,y):
        self.delta = -( y - self.Z ) * self.grad_Z
        return self.delta
        
class NeuroNetwork(Classifier):
    def __init__(self,topology,K,seuil,eta=0.2,iterNum=200):
        self.err = []
        self.topology    = topology
        self.hiddenLayers= np.array([ HiddenLayer(inputSize,layerSize,g,grad_g,eta) for g,grad_g,inputSize,layerSize in self.topology])
        self.outputLayer = OutputLayer(layerSize,K,g,grad_g,eta)
        self.layers      = np.append(self.hiddenLayers,self.outputLayer)
        self.eta         = eta
        self.iterNum     = iterNum
        self.seuil       = seuil
    def fit(self,x,y):
        self.x   = np.array( [ np.append(line,1) for line in x ] )
        self.y   = y
        self.optimize()
    def forward(self,X):
        Z = X
        for l in self.layers:
            Z = l.forward(Z)
        return Z
    def backProp(self,x,y):
        # errors
        deltaWeiSum = y
        for i in range(len(self.layers)-1,-1,-1):
            delta       = self.layers[i].calculateDelta(deltaWeiSum)
            wm          = self.layers[i].getWeightMatrix()
            deltaWeiSum = np.sum(delta * wm.T,axis=1)

        # update layer weights
        for i in range(len(self.layers)):
            self.layers[i].updateWeights()

    def optimize(self):
        for i in range(len(self.x)):
            self.forward(self.x[i])
            self.err.append((self.y[i] - self.outputLayer.Z)**2/2.)
            self.backProp(self.x[i],self.y[i])
        
    def predict(self,X):
        resultats = []
        for x in X:
            x = np.append(x,1)
            self.forward(x)
            if self.outputLayer.Z >= self.seuil:
                resultats.append(1.)
            else:
                resultats.append(-1.)
        return np.array(resultats)
        
trainX,trainY = gen_arti(data_type=3)
testX ,testY  = gen_arti(data_type=3)
# np.save("trainX",trainX[:10])
# np.save("trainY",trainY[:10])
# trainX = np.load("trainX.npy")
# trainY = np.load("trainY.npy")
# trainY = np.array(zip(trainY,trainY))
eta = 0.2
ann = NeuroNetwork([[sigmoid,grad_sigmoid,3,10],[sigmoid,grad_sigmoid,10,6]],1,0.5,eta)
# ann = NeuroNetwork([[sigmoid,grad_sigmoid,3,15]],1,0.5,eta)
# ann = NeuroNetwork([[tangH,grad_tangH,3,5]],1,eta)

# print ann.layers[-1].getWeightMatrix()

while ann.score(testX,testY) < 0.9:
    ann.fit(trainX,trainY)
    print ann.score(testX,testY)

# resultats = ann.predict([ np.append(x,1) for x in testX ])
'''
fig = plt.figure()
plt.plot(ann.err)
plt.show()
plt.close(fig)
# print ann.score(testX,testY)
plot_frontiere(trainX,ann.predict)
plot_data(trainX,trainY)
plt.show()
'''
