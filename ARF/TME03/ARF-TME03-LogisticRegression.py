#!/bin/python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from arftools import *

class LogisticRegression(Classifier):
    def __init__(self,eps=1e-4,max_iter=5000,delta=1e-6,dim=2):
        self.dim = dim
        GradientDescent.__init__(self,self,eps,max_iter,delta)
    def fit(self,x,y):
        self.data = data
        self.y    = y
        self.optimize()
    def f(self,w):
        return 
    def grad_f(self,w):
        return 
    def predict(self,x):
        X = to_array(testX)
        return self.f()


