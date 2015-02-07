# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # TME1
# 
# _Objectifs_ :
# 
# * prendre en main la syntaxe objet python (et les modules)
# 
# * expérimenter les arbres de décisions
# 
# * découverte des effets de sur/sous-apprentissage
# 
# 
# ## Arbre de décisions et objets

# <markdowncell>

# Le code suivant est une implémentation objet des arbres de décision en python, pour variable continue.
# Il y a 3 passages à compléter (une ligne à chaque fois) :
# 
# * dans la classe mère Classifier, la fonction score(x,y), qui doit donner la précision sur l'ensemble (x,y) en utilisant la fonction predict
# 
# * dans la classe DecisionTree, dans fit et predict
# 
# Lisez bien le code de facon a comprendre grossierement comment fonctionne les objets en python. Completez ce qui est necessaire.

# <codecell>

import numpy as np
from collections import Counter
import pydot  #pour l'affichage graphique d'arbres

###############################
# Fonctions auxiliaires     
###############################

def p_log_p(counts):
    """ fonction pour calculer \sum p_i log(p_i) """
    return np.nan_to_num(np.sum(counts*np.log2(counts)))

def entropy(y):
    """ calcul de l'entropie d'un ensemble, attention c'est lent! """
    ylen = float(y.size)
    if ylen <= 1:
        return 0
    counts = np.array(Counter(y).values())/ylen
    return -p_log_p(counts)

###############################
# Classes
###############################


class Classifier(object):
    """ Classe generique d'un classifieur
        Dispose de 3 méthodes : 
            fit pour apprendre
            predict pour predire
            score pour evaluer la precision
    """
    
    def fit(self,x,y):
        raise NotImplementedError("fit non  implemente")
    def predict(self,x):
        raise NotImplementedError("predict non implemente")
    def score(self,x,y):
        return (self.predict(x) == y).mean()

class Split(object):
    """ Permet de coder un split pour une variable continue
        Contient :
            * le numero de la variable ou se fait le split
            * le seuil du split
            * le gain d'information du split
         predict(x) renvoie -1 si x_i<=seuil, +1 sinon
         best_gain(x,y) calcul le meilleur seuil pour la colonne x (1-dimension) et les labels y
         find_best_split(x,y) calcul le meilleur split pour les données x (n-dimensions) et les labels y
     """
    def __init__(self,idvar=None,threshold=None,gain=None):
        self.idvar=idvar
        self.threshold=threshold
        self.gain=gain

    def predict(self,x):
        if len(x.shape)==1:
            x=x.reshape((1,x.shape[0]))
        return [-1 if x[i,self.idvar]<=self.threshold else 1 for i in range(x.shape[0])]

    @staticmethod
    def best_gain(x,y):
        ylen = float(y.size)
        idx_sorted = np.argsort(x)
        h=entropy(y)
        xlast=x[idx_sorted[0]]
        split_val=x[idx_sorted[0]]
        hmin = h
        for i in range(y.size):
            if x[idx_sorted[i]]!=xlast:
                htmp = i/ylen*entropy(y[idx_sorted[:i]])+(ylen-i)/ylen*entropy(y[idx_sorted[i:]])
                if htmp<hmin:
                    hmin=htmp
                    split_val=(xlast+x[idx_sorted[i]])/2.
            xlast=x[idx_sorted[i]]
        return (h-hmin/ylen),split_val

    @staticmethod
    def find_best_split(x,y):
        if len(x.shape)==1:
            x = x.reshape((1,x.shape[0]))
        hlist = [[Split.best_gain(x[:,i],y),i] for i in range(x.shape[1])]
        (h,threshold),idx= max(hlist)
        return Split(idx,threshold,h)

    def __str__(self):
        return "var %s, thresh %f (gain %f)" %(self.idvar,self.threshold, self.gain)

class Node(Classifier):
    """ Noeud d'un arbre
        split contient le split du noeud
        parent,left, right les noeuds parents, gauche et droit du noeud (parent = None si racine, left=right=None si feuille)
        leaf indique si le noeud est une feuille
        depth la profondeur du noeud
        label le label majoritaire du noeud
        info est un dictionnaire pour stocker des informations supplemantaires
    """
    def __init__(self,split=None,parent=None,left=None,right=None,leaf=True,depth=-1,label=None,**kwargs):
        self.split=split
        self.parent=None
        self.left=None
        self.right=None
        self.info=dict(kwargs)
        self.leaf=leaf
        self.label=label
        self.depth=depth

    def predict(self,x):
        if len(x.shape)==1:
            x=x.reshape((1,x.shape[0]))
        if self.leaf:
            return [self.label]*x.shape[0]
        return [self.left.predict(x[i,:])[0] if res<0 else self.right.predict(x[i,:])[0] for i,res in enumerate(self.split.predict(x))]


    def fit(self,x,y):
        counts=Counter(y)
        self.split=Split.find_best_split(x,y)
        self.label = counts.most_common()[0][0]

    def __str__(self):
        if self.leaf:
            return "Leaf : %s" % (self.label,)
        return "Node : %s (%s)" % (self.split,self.info)

class DecisionTree(Classifier):
    """ Arbre de decision
    max_depth indique la profondeur max de l'arbre
    min_samples_split le nombre minimal d'exemples pour continuer l'apprentissage
    to_dot permet de convertir en dot l'arbre (affichage graphique)
    to_pdf d'enregistrer l'arbre dans un fichier pdf
    """
    
    def __init__(self,max_depth=None,min_samples_split=2):
        self.max_depth=max_depth
        self.min_samples_split=min_samples_split

    def fit(self,x,y):
        """ apprentissage de l'arbre de maniere iterative
        on apprend un noeud, puis on cree les deux enfants de ce noeud, que l'on ajoute a la pile des noeuds
        a traiter par la suite (nodes_to_treat), ainsi que les index des exemples associes (dic_idx)
        """
        
        self.root=Node(depth=0)
        nodes_to_treat = [self.root]
        dic_idx=dict({self.root : range(len(y))})
        while len(nodes_to_treat)>0:
            # recuperation du noeud courant
            curnode = nodes_to_treat.pop()
            #recuperation de la liste des indices des exemples associes, x[idx_train,:] contient l'ensemble des 
            #exemples a traiter
            idx_train = dic_idx.pop(curnode)
            # infos complementaires sur le nombre d'exemples en apprentissage par label
            for lab,clab in Counter(y[idx_train]).items():
                curnode.info[lab]=clab
            
            # A COMPLETER #
            #trouve le meilleur split pour ce noeud
            curnode.fit(x[idx_train],y[idx_train])
            
            # recupere les predictions pour partager entre fils droit et gauche les exemples 
            pred =  curnode.split.predict(x[idx_train,:])
            l_idx = [ idx_train[i] for i in range(len(idx_train)) if pred[i]<0 ]
            r_idx = list(set(idx_train).difference(l_idx))
            
            #Condition d'arrets
            if entropy(y[idx_train])==0 or curnode.depth >= self.max_depth or \
                    len(l_idx) < self.min_samples_split or len(r_idx) < self.min_samples_split:
                curnode.leaf=True
                continue
            #Creation des deux enfants
            curnode.left = Node(parent=curnode,depth=curnode.depth+1)
            curnode.right = Node(parent=curnode,depth=curnode.depth+1)
            curnode.leaf=False
            #On enregistre les indices correspondant aux deux noeuds
            dic_idx[curnode.left]=l_idx
            dic_idx[curnode.right]=r_idx
            #On ajoute les deux enfants a la liste des noeuds a traiter
            nodes_to_treat = [curnode.left,curnode.right]+nodes_to_treat
        
    def predict(self,x):
        decision = self.root
        while decision.leaf == False:
            sum(decision.predict(x)) 
            
            
    
        
    def __str__(self):
        s=""
        nodes=[self.root]
        while len(nodes)>0:
            curnode=nodes.pop()
            if not curnode.leaf:
                s+= "\t"*curnode.depth + "var %d :  <=|> %f \n"  %(curnode.split.idvar,curnode.split.threshold)
                nodes+=[curnode.left,curnode.right]
            else:
                s+= "\t"*curnode.depth + "class : %s\n" %(curnode.label,)
        return s

    def to_dot(self,dic_var=None):
        s="digraph Tree {"
        cpt=0
        nodes = [(self.root,cpt)]
        while len(nodes)>0:
            curnode,idx = nodes.pop()
            labinfo = ",".join(["%s: %s" % (lab,slab) for lab,slab in curnode.info.items()])
            if not curnode.leaf:
                s+="%d [label=\"%s <= %f\n IG=%f\n " %(idx,curnode.split.idvar \
                    if not dic_var else dic_var[curnode.split.idvar],curnode.split.threshold,curnode.split.gain)
                s+= " %s \n \",shape=\"box\" ];\n"  % (labinfo,)
                lidx = cpt +1
                ridx = cpt +2
                s+= "%d -> %d; %d -> %d;\n" % (idx,lidx,idx,ridx)
                cpt+=2
                nodes += [(curnode.left,lidx),(curnode.right,ridx)]
            else:
                s+= "%d [label=\"label=%s\n %s \"];\n" %(idx,curnode.label,labinfo)
        return s+"}"
    def to_pdf(self,filename,dic_var=None):
        pydot.graph_from_dot_data(self.to_dot(dic_var)).write_pdf(filename)

# <markdowncell>

# ## Expérimentations sur USPS
# 
# Tester l'algorithme sur les données du [TME3 de MAPSI](http://webia.lip6.fr/~mapsi/pmwiki.php?n=Cours.TME3). Servez vous soit du code du TME3, soit du suivant pour lire le fichier de données. Attention ! L'implémentation est lourde, tester l'algorithme sur un arbre de petite profondeur. 
# 
# Visualisez sur quelques exemples les arbres. Observez l'erreur sur l'ensemble d'apprentissage. Comment se comporte-t-elle en fonction des deux paramètres ? Est-elle une bonne prédiction de l'erreur de votre modèle ? Comment de manière simple obtenir une meilleure prédiction (ensemble de test) ?

# <codecell>

def load_usps(filename):
    with open(filename,"r") as f:
        f.readline()
        data =[ [float(x) for x in l.split()] for l in f if len(l.split())>2]
    tmp = np.array(data)
    return tmp[:,1:],tmp[:,0].astype(int)

# <markdowncell>

# Choisissez dans la suite quelques classes parmi celles disponibles. Partager votre ensemble en deux sous-ensembles, un d'apprentissage qui vous servira à apprendre votre modèle, l'autre de test qui vous servira à évaluer l'erreur.
# 
# + Tracez en fonction de la profondeur l'erreur en apprentissage et en test
# 
# + Faites varier la taille de vos deux ensembles. Que remarquez-vous ?
# 
# 
# ## Expérimentations sur jeu de données artificielles 
# 
# Que font les fonctions randCheckers et bigauss ?
# Recommencez vos expériences sur ces 2 jeux de données. Que remarquez-vous ?
# 

# <codecell>

def randCheckers(n1,n2,epsilon=0.1):
    nbp=int(numpy.floor(n1/8))
    nbn=int(numpy.floor(n2/8))
    xapp=np.reshape(random.rand((nbp+nbn)*16),[(nbp+nbn)*8,2])
    yapp=[1]*((nbp+nbn)*8)
    idx=0
    for i in range(-2,2):
        for j in range(-2,2):
            if (((i+j) % 2)==0):
                nb=nbp
            else:
                nb=nbn
                yapp[idx:(idx+nb)]=[-1]*nb
            xapp[idx:(idx+nb),0]=np.random.rand(nb)+i+epsilon*random.randn(nb)
            xapp[idx:(idx+nb),1]=np.random.rand(nb)+j+epsilon*random.randn(nb)
            idx=idx+nb
    ind=range((nbp+nbn)*8)
    random.shuffle(ind)
    return xapp[ind,:],yapp[ind]
  
def bigauss(n,mu1=[1,1],mu2=[-1,-1],sigma1=[0.2,0.2],sigma2=[0.5,0.5]):
    x=np.vstack((np.random.multivariate_normal(mu1,np.diag(sigma1),n),np.random.multivariate_normal(mu2,np.diag(sigma2),n)))
    y=np.vstack((np.ones((n,1)),-np.ones((n,1))))
    ind=np.random.permutation(range(2*n))
    return x[ind,:],y[ind,:]



# <codecell>


