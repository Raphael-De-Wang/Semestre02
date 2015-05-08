# -*- coding: utf-8 -*-
#!env python

from ANN import *
'''
print 'test input : ', np.arange(-5,6)/2.
print 'test toInt : ', map(toInt,np.arange(-5,6)/2.)
print 'test a_unif: ', np.array(map(a_unif,np.arange(41),np.ones(41)*41))
print 'test b_unif: ', np.array(map(b_unif,np.arange(11),np.ones(11)*11))
print 'test c_unif: ', np.array(map(c_unif,np.arange(11),np.ones(11)*11))
'''

trainX1,trainY1 = data_uniform()
trainX2,trainY2 = data_random(data_size=10000-len(trainY1))
trainX = np.vstack([trainX1,trainX2])
trainY = np.vstack([trainY1,trainY2])
# trainX ,trainY  = data_random(data_size=1000)
testX ,testY  = data_random(data_size=10)

eta = 0.01
ann = NeuroNetwork([[sigmoid,grad_sigmoid,64,24]],91,0.5,eta)
for i in range(4):
    ann.fit(trainX,trainY)
    print ann.predict(testX) - 45
    print np.argmax(testY,axis=1) - 45
        

