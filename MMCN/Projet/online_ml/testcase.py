# -*- coding: utf-8 -*-
#!env python

from ANN import *

'''
print 'test a_unif: ', np.array(map(a_unif,np.arange(41),np.ones(41)*41))
print 'test input : ', np.arange(-5,6)/2.
print 'test toInt : ', map(toInt,np.arange(-5,6)/2.)
print 'test b_unif: ', np.array(map(b_unif,np.arange(11),np.ones(11)*11))
print 'test c_unif: ', np.array(map(c_unif,np.arange(11),np.ones(11)*11))
'''

N1  = 41
N2  = 11
N3  = 11
N   = 10000

trainX1,trainY1 = data_uniform()
trainX2,trainY2 = data_random(data_size=N-len(trainY1))
trainX = np.vstack([trainX1,trainX2])
trainY = np.vstack([trainY1,trainY2])

indice = range(N)
# np.random.shuffle(indice)
# trainX ,trainY  = data_random(data_size=100)
testX ,testY  = data_random(data_size=10)

ann = NeuroNetwork([[sigmoid,grad_sigmoid,64,21]],91,eta=0.1)
for i in range(1):
    ann.fit(trainX[indice],trainY[indice])
    # print np.argmax(testY,axis=1) - 45
    yk, ykCible = ann.predict(testX)
    print "ykCible : ",ykCible
    print "yk      : ",yk
        

