# -*- coding: utf-8 -*-
#!env python

from   ANN               import *
import matplotlib.pyplot as     plt
from operator import itemgetter

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

# ann = NeuroNetwork([[sigmoid,grad_sigmoid,64,21]],91,eta=0.01)
# ann.fit(trainX,trainY)

# yk, ykCible = ann.predict(testX)
# print "ykCible   : ",ykCible
# print "yk        : ",yk
'''
yk, ykCible = ann.predict(testX)
print "ykCible   : ",ykCible
print "yk        : ",yk
'''

# 2. Représenter graphiquement les courbes d'accord
def tuningCurves(fname=None):
    fig = plt.figure(figsize=(20, 15))
    plt.suptitle('Tuning Curves', fontsize=35)
    x = 25
    r = 25
    plt.subplot(331,axisbg="#fdf6e3")
    rObj = [ posAct(x,i,N1) for i in range(N1) ]
    pos = np.argmax(rObj)
    plt.plot(rObj,color="#2aa198",linewidth="2")
    plt.plot(pos*np.ones(11),np.arange(0,1.05,0.1),'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity x=%d r=%d"%(x,r))
    
    plt.subplot(332,axisbg="#fdf6e3")
    rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
    plt.plot(range(N1,N1+N2),rPlu,color="#2aa198",linewidth="2")
    plt.ylim(ymax=1)
    
    plt.subplot(333,axisbg="#fdf6e3")
    rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
    plt.plot(range(N1+N2,N1+N2+N3),rMin,color="#2aa198",linewidth="2")
    plt.ylim(ymax=1)
    
    x = 0
    r = 0
    plt.subplot(334,axisbg="#fdf6e3")
    rObj = [ posAct(x,i,N1) for i in range(N1) ]
    pos = np.argmax(rObj)
    plt.plot(rObj,color="#2aa198",linewidth="2")
    plt.plot(pos*np.ones(11),np.arange(0,1.05,0.1),'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity x=%d r=%d"%(x,r))
    
    plt.subplot(335,axisbg="#fdf6e3")
    rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
    plt.plot(range(N1,N1+N2),rPlu,color="#2aa198",linewidth="2")
    plt.ylim(ymax=1)

    plt.subplot(336,axisbg="#fdf6e3")
    rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
    plt.plot(range(N1+N2,N1+N2+N3),rMin,color="#2aa198",linewidth="2")
    plt.ylim(ymax=1)

    x = -25
    r = -25
    plt.subplot(337,axisbg="#fdf6e3")
    rObj = [ posAct(x,i,N1) for i in range(N1) ]
    pos = np.argmax(rObj)
    plt.plot(rObj,color="#2aa198",linewidth="2")
    plt.plot(pos*np.ones(11),np.arange(0,1.05,0.1),'--',color="#cb4b16",linewidth="2")
    plt.xlabel("Position Neuron Index")
    plt.ylabel("Activity x=%d r=%d"%(x,r))
    
    plt.subplot(338,axisbg="#fdf6e3")
    rPlu = [ rPlAct(r,i,N2) for i in range(N2) ]
    plt.plot(range(N1,N1+N2),rPlu,color="#2aa198",linewidth="2")
    plt.xlabel("Direction Left to Right Neuron Index")
    plt.ylim(ymax=1)
    
    plt.subplot(339,axisbg="#fdf6e3")
    rMin = [ rMiAct(r,i,N3) for i in range(N3) ]
    plt.plot(range(N1+N2,N1+N2+N3),rMin,color="#2aa198",linewidth="2")
    plt.xlabel("Direction Right to Left Neuron Index")
    plt.ylim(ymax=1)
    if fname==None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)

# 3. Représenter graphiquement la sortie du réseau
def graph_sorti(ann,minErrs,hlNum=21,fname=None):
    fig = plt.figure(figsize=(20, 15))
    plt.suptitle('Sortie du Reseau, Hidden Layer Size : %d'%hlNum, fontsize=35)
    plt.subplot(311,axisbg="#fdf6e3")
    i = minErrs[0]
    Y = ann.get_output(testX[i])
    y = np.arange(np.max(Y),-0.2,-0.01)
    x = np.argmax(Y)*np.ones(len(y))
    plt.plot(Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="2")
    plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity")
    plt.legend()
    
    plt.subplot(312,axisbg="#fdf6e3")
    i = minErrs[1]
    Y = ann.get_output(testX[i])
    y = np.arange(np.max(Y),-0.2,-0.01)
    x = np.argmax(Y)*np.ones(len(y))
    plt.plot(Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="2")
    plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity")
    plt.legend()
    
    plt.subplot(313,axisbg="#fdf6e3")
    i = minErrs[2]
    Y = ann.get_output(testX[i])
    y = np.arange(np.max(Y),-0.2,-0.01)
    x = np.argmax(Y)*np.ones(len(y))
    plt.plot(Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="2")
    plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.xlabel("Neuron Index")
    plt.ylabel("Activity")
    plt.legend()
    if fname==None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)

# 3. Représenter graphiquement la couche cachée du réseau
def graph_hidden(ann,minErrs,hlNum=21,fname=None):
    fig = plt.figure(figsize=(20, 15))
    plt.suptitle('Sortie de la Couche Cachee du Reseau, Hidden Layer Size : %d'%hlNum, fontsize=35)
    plt.subplot(311,axisbg="#fdf6e3")
    i = minErrs[0]
    Y = ann.get_hiddenZ(testX[i])
    X = range(21)
    # y = np.arange(np.max(Y),-0.2,-0.01)
    # x = np.argmax(Y)*np.ones(len(y))
    plt.bar(X,Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="1")
    # plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity")
    plt.legend()
    
    plt.subplot(312,axisbg="#fdf6e3")
    i = minErrs[10]
    Y = ann.get_hiddenZ(testX[i])
    # y = np.arange(np.max(Y),-0.2,-0.01)
    # x = np.argmax(Y)*np.ones(len(y))
    plt.bar(X,Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="1")
    # plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.ylabel("Activity")
    plt.legend()
    
    plt.subplot(313,axisbg="#fdf6e3")
    i = minErrs[200]
    Y = ann.get_hiddenZ(testX[i])
    # y = np.arange(np.max(Y),-0.2,-0.01)
    # x = np.argmax(Y)*np.ones(len(y))
    plt.bar(X,Y, label="x=%d, r=%d"%(testX[i][0],testX[i][-1]),color="#2aa198",linewidth="1")
    # plt.plot(x,y,'--',color="#cb4b16",linewidth="2")
    plt.xlabel("Neuron Index")
    plt.ylabel("Activity")
    plt.legend()
    if fname==None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)

# 4. Estimer la performance du réseau
def performance(err_list,hlNum=21,fname=None):
    moy   = np.mean(err_list)
    moyAbs= np.mean(abs(err_list))
    std   = np.std(err_list)
    print "Erreur Moyenne : ", moy
    print "Erreur Moyenne Abs: ", moyAbs
    print "L'écart-type de l'erreur : ", std
    fig = plt.figure(figsize=(20, 15))
    ax = fig.add_subplot(111,axisbg="#fdf6e3")
    plt.suptitle('Performance du Reseau, Hidden Layer Size : %d'%hlNum, fontsize=35)
    plt.hist(err_list.tolist(),color="#859900")
    plt.legend(["err_moy=%f\nerr_moy_abs=%f\necart-type=%f"%(moy,moyAbs,std)],prop={'size':16})
    plt.ylabel("Occurences")
    plt.xlabel("Erreur")
    if fname==None:
        plt.show()
    else:
        plt.savefig(fname)
    plt.close(fig)
    return moy,moyAbs,std

tuningCurves("TuningCurves.png")

# training ANN et estimer les erreurs
y_vrai = lambda x : x[0] + x[-1]
erreur = lambda yv, ye : yv - ye

ann = None

for hlNum in [21]:
    annList    = []
    moyList    = []
    moyAbsList = []
    stdList    = []
    for i in range(1):
        elist  = []
        minErrs= []
        ann = NeuroNetwork([[sigmoid,grad_sigmoid,64,hlNum]],91,eta=0.01)
        ann.fit(trainX,trainY)
        for x in testX:
            yv = y_vrai(x)
            err= erreur(yv,ann.y_esti(x))
            elist.append(err)
        elist = np.array(elist)
        minErrs = np.array(sorted(zip(np.arange(len(elist)),np.abs(elist)),key=itemgetter(1)))[:,0]
        graph_sorti(ann,minErrs,hlNum,"LaSortieDuReseau[%d][hlNum%d].png"%(i,hlNum))
        graph_hidden(ann,minErrs,hlNum,"LaSortieCoucheCacheeDuReseau[%d][hlNum%d].png"%(i,hlNum))
        moy,moyAbs,std = performance(elist,hlNum,"LaPerformanceDuReseau[%d][hlNum%d].png"%(i,hlNum))
        moyList.append(moy)
        moyAbsList.append(moyAbs)
        stdList.append(std)
    # np.savetxt("err_moy_std[hlNum%d].txt"%hlNum,np.array(zip(moyList,moyAbsList,stdList)))


    
