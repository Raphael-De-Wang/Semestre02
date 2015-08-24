#!env python

import csv
import sys
import pylab
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plotgCAIsDomain(X,Y,xDiv=1.,yDiv=20.,figName=None):
    fig = pylab.figure(figsize=(320,3))
    # pylab.hist2d(X,Y,bins=[len(X)/xDiv,yDiv], norm=mpl.colors.LogNorm(), cmap=mpl.cm.jet)
    pylab.hist2d(X,Y,bins=[len(np.unique(X))/xDiv,yDiv],norm=mpl.colors.LogNorm(), cmap=mpl.cm.jet)
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.title("Desity Plot - Domains Grouped by [%d], gCAIs Divided by [%d], Top [%d] Domains (PS : values in brackets are manipulable)"%(xDiv,yDiv,len(np.unique(Y))))
    fig.subplots_adjust(bottom=0.25)
    cbar = pylab.colorbar()
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close(fig)

def plotgCAIsDomain2(X,Y,figName=None,xticks=[]):
    T = np.arctan(Y) #,np.ones(len(X)))
    fig = pylab.figure(figsize=(30,4))
    xUniLen = len(np.unique(X))
    pylab.scatter(X, Y, s=5, c=T, alpha=.4, marker='x')
    pylab.xlabel("Domain Sorted By Translation Level")
    pylab.ylabel("gCAIs")
    pylab.ylim(0,1.02), pylab.yticks(pylab.arange(0,1.1,0.1))
    pylab.xlim(0,xUniLen), pylab.xticks(xticks)
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen)*0.025,'k-')
    pylab.plot(np.arange(xUniLen),np.ones(xUniLen),'k-')
    pylab.title("gCAIs values distributions ordered by expression level, the colors correspond to gCAI value, each colume represent a domain")
    if figName:
        pylab.savefig(figName)
    else:
        pylab.show()
    pylab.close()
    
def load_domain_gcai_abundance():
    x = 0
    rsv_dname = ""
    X = []
    Y = []
    dlist = []
    with open("domain_gcai_abundance.csv") as csv_handle :
        reader = csv.reader(csv_handle, delimiter='\t')
        for dname, gcai, abundance in reader :
            if cmp(rsv_dname, dname) <> 0:
                x += 1
                rsv_dname = dname
                dlist.append(dname)
            X.append(x)
            Y.append(float(gcai))
    return (X,Y,dlist)
    
(X,Y,dlist) = load_domain_gcai_abundance()
plotgCAIsDomain(X,Y,figName="gCAI_order_by_abundance")
plotgCAIsDomain2(X,Y,figName="gCAI_order_by_abundance2")

