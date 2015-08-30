#!env python

import matplotlib.cm as cm
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,40))

for p in range(11):
    ax = plt.subplot(11,3,p*3+1)
    ax.set_title("[x, y, p=%.1f]"%(p/10.))
    X = [ x/10. for y in range(11) for x in range(11) ]
    Y = [ y/10. for y in range(11) for x in range(11) ]
    t = [ (x/10.,y/10.,p/10.,1) for y in range(11) for x in range(11) ]
    plt.scatter(X,Y,c=t,s=300)

for p in range(11):    
    ax = plt.subplot(11,3,p*3+2)
    ax.set_title("[x, p=%.1f, y]"%(p/10.))
    X = [ x/10. for y in range(11) for x in range(11) ]
    Y = [ y/10. for y in range(11) for x in range(11) ]
    t = [ (x/10.,p/10.,y/10.,1) for y in range(11) for x in range(11) ]
    plt.scatter(X,Y,c=t,s=300)

for p in range(11):    
    ax = plt.subplot(11,3,(p+1)*3)
    ax.set_title("[p=%.1f, x, y]"%(p/10.))
    X = [ x/10. for y in range(11) for x in range(11) ]
    Y = [ y/10. for y in range(11) for x in range(11) ]
    t = [ (p/10.,x/10.,y/10.,1) for y in range(11) for x in range(11) ]
    plt.scatter(X,Y,c=t,s=300)

# plt.show()
plt.savefig("/tmp/colormap.pdf")
