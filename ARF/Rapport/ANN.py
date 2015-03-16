# -*- coding: utf-8 -*-
#!env python

import numpy as np
import matplotlib.pyplot as plt

from arftools import *

def sigmoid(Neuro):
    pass
    
class Neuron(OptimFunc):
    def x_random(self):
        return np.ones(self.dim)/self.dim
    def a(self,x):
        return sum(self.weight*x)
    def f(self,x):
        return self._f(self.a(x))
    def grad_f(self,x):
        return self._grad_f(self.a(x))
        
class Layer(object):
    def __init__(self,inputSize,g,layerSize,layerName):
        self.g         = g
        self.inputSize = inputSize
        self.layerSize = layerSize
        self.layerName = layerName
        self.layer     = np.array([ Neuron(inputSize,g) for i in num ])
        self.bias      = Neuron(inputSize,lambda:x return 1,lambda:x return 1)
        self.layer     = np.append(self.layer,bias)
    def forword(self,prevZ):
        return np.array([ n.forward(prevZ) for n in self.layer ])
        
class NeuroNetwork(Classifier,OptimFunc,GradientDescent):
    def __init__(self,topology,eps=1e-4,max_iter=5000,delta=1e-6):
        self.topology = topology
        GradientDescent.__init__(self,self,eps,max_iter,delta)
    def fit(self,x,y):
        self.dim      = len(data[0])
        self.topology[0,1] = self.dim
        self.layers   = np.array([ Layer(inputSize,g,layerSize,"Hidden Layer") for g,inputSize,layerSize in self.topology])
        self.data     = data
        self.y        = y
        self.optimize()
    def forward(self):
        Z = []
        for prevZ in self.x:
            for l in self.layers:
                prevZ = l.forword(prevZ)
            Z.append(prevZ)
        return np.array(Z)
        
    def backProp(self):
        pass
    def predict(self,x):
        pass

