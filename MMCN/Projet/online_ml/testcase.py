# -*- coding: utf-8 -*-
#!env python

from   ANN               import *
import matplotlib.pyplot as     plt

'''
print 'test a_unif: ', np.array(map(a_unif,np.arange(41),np.ones(41)*41))
print 'test input : ', np.arange(-5,6)/2.
print 'test toInt : ', map(toInt,np.arange(-5,6)/2.)
print 'test b_unif: ', np.array(map(b_unif,np.arange(11),np.ones(11)*11))
print 'test c_unif: ', np.array(map(c_unif,np.arange(11),np.ones(11)*11))
'''

# 1. Ecrire un programme

N1  = 41
N2  = 11
N3  = 11
N   = 10000

trainX1,trainY1 = data_uniform()
trainX2,trainY2 = data_random(data_size=N-len(trainY1))
trainX = np.vstack([trainX1,trainX2])
trainY = np.vstack([trainY1,trainY2])

# trainX,trainY = data_random(data_size=100)
testX ,testY  = data_random(data_size=1000)

ann = NeuroNetwork([[sigmoid,grad_sigmoid,64,21]],91,eta=0.01)
ann.fit(trainX,trainY)
'''
yk, ykCible = ann.predict(testX)
print "ykCible   : ",ykCible
print "yk        : ",yk
'''

# 2. Représenter graphiquement les courbes d'accord
fig = plt.figure()
x = 25
r = 25
plt.subplot(331)
rObj = [ posAct(x,i,N1) for i in range(N1) ]
plt.plot(rObj)
plt.ylabel("Activity x=%d r=%d"%(x,r))

plt.subplot(332)
rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
plt.plot(range(N1,N1+N2),rPlu)
plt.ylim(ymax=1)

plt.subplot(333)
rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
plt.plot(range(N1+N2,N1+N2+N3),rMin)
plt.ylim(ymax=1)

x = 0
r = 0
plt.subplot(334)
rObj = [ posAct(x,i,N1) for i in range(N1) ]
plt.plot(rObj)
plt.ylabel("Activity x=%d r=%d"%(x,r))

plt.subplot(335)
rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
plt.plot(range(N1,N1+N2),rPlu)
plt.ylim(ymax=1)

plt.subplot(336)
rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
plt.plot(range(N1+N2,N1+N2+N3),rMin)
plt.ylim(ymax=1)

x = -25
r = -25
plt.subplot(337)
rObj = [ posAct(x,i,N1) for i in range(N1) ]
plt.plot(rObj)
plt.xlabel("Position Neuron Index")
plt.ylabel("Activity x=%d r=%d"%(x,r))

plt.subplot(338)
rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
plt.plot(range(N1,N1+N2),rPlu)
plt.xlabel("Direction Left to Right Neuron Index")
plt.ylim(ymax=1)

plt.subplot(339)
rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
plt.plot(range(N1+N2,N1+N2+N3),rMin)
plt.xlabel("Direction Right to Left Neuron Index")
plt.ylim(ymax=1)

plt.show()
plt.close(fig)

# 3. Représenter graphiquement la sortie du réseau
fig = plt.figure()
plt.subplot(311)
i = 1
plt.plot(ann.get_output(testX[i]), label="x=%d, r=%d"%(testX[i][0],testX[i][-1]))
plt.ylabel("Activity")
plt.legend()

plt.subplot(312)
i = 2
plt.plot(ann.get_output(testX[i]), label="x=%d, r=%d"%(testX[i][0],testX[i][-1]))
plt.ylabel("Activity")
plt.legend()

plt.subplot(313)
i = 3
plt.plot(ann.get_output(testX[i]), label="x=%d, r=%d"%(testX[i][0],testX[i][-1]))
plt.xlabel("Neuron Index")
plt.ylabel("Activity")
plt.legend()

plt.show()
plt.close(fig)

# 4. Estimer la performance du réseau
y_vrai = lambda x : x[0] + x[-1]
erreur = lambda yv, ye : yv - ye
elist  = []
for x in testX:
    yv = y_vrai(x)
    err= erreur(yv,ann.y_esti(x))
    elist.append(err)
fig = plt.figure()
plt.hist(elist)
plt.show()
plt.close(fig)
