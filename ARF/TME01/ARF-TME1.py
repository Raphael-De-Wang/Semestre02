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
import matplotlib.pyplot as plt
###############################
# Fonctions auxiliaires
###############################

eps = np.finfo(float).eps

def p_log_p(counts):
    """ fonction pour calculer \sum p_i log(p_i) """
    return np.sum(counts*np.log2(counts+eps))

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
       """ A COMPLETER """
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
        return np.array([-1 if x[i,self.idvar]<=self.threshold else 1 for i in range(x.shape[0])])

    @staticmethod
    def best_gain(x,y):
        ylen = float(y.size)
        idx_sorted = np.argsort(x)
        h=entropy(y)
        counts = Counter(y)
        dic_labs=dict(zip(counts.keys(),range(len(counts))))
        rcount = np.array(counts.values())
        lcount = np.zeros(len(counts))
        lc=0.
        rc=ylen
        xlast=x[idx_sorted[0]]
        split_val=x[idx_sorted[0]]
        hmin = h
        for i in range(y.size):
            if x[idx_sorted[i]]!=xlast:
                #htmp = i/ylen*entropy(y[idx_sorted[:i]])+(ylen-i)/ylen*entropy(y[idx_sorted[i:]])
                htmp = i /ylen*(-p_log_p(lcount/lc))+(ylen-i)/ylen*(-p_log_p(rcount/rc))
                if htmp<hmin:
                    hmin=htmp
                    split_val=(xlast+x[idx_sorted[i]])/2.
            rcount[dic_labs[y[idx_sorted[i]]]]-=1
            lcount[dic_labs[y[idx_sorted[i]]]]+=1
            lc+=1
            rc+=1
            xlast=x[idx_sorted[i]]
        return (h-hmin),split_val

    @staticmethod
    def find_best_split(x,y):
        if len(x.shape)==1:
            x = x.reshape((1,x.shape[0]))
        hlist = [[Split.best_gain(x[:,i],y),i] for i in range(x.shape[1])]
        (h,threshold),idx= max(hlist)
        return Split(idx,threshold,h)

    def __str__(self):
        return "var %s, thresh %f (gain %f)" %(self.idvar,self.threshold, self.gain)
    def __repr__(self):
        return self.__str__()

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
            return np.array([self.label]*x.shape[0])
        return np.array([self.left.predict(x[i,:])[0] if res<0 else self.right.predict(x[i,:])[0] for i,res in enumerate(self.split.predict(x))])

    def fit(self,x,y):
        counts=Counter(y)
        self.split=Split.find_best_split(x,y)
        self.label = counts.most_common()[0][0]

    def __str__(self):
        if self.leaf:
            return "Leaf : %s" % (self.label,)
        return "Node : %s (%s)" % (self.split,self.info)
    def __repr__(self):
        return self.__str__()

class DecisionTree(Classifier):
    """ Arbre de decision
    max_depth indique la profondeur max de l'arbre
    min_samples_split le nombre minimal d'exemples pour continuer l'apprentissage
    to_dot permet de convertir en dot l'arbre (affichage graphique)
    to_pdf d'enregistrer l'arbre dans un fichier pdf
    """

    def __init__(self,max_depth=5,min_samples_split=2):
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
        # A COMPLETER
        return self.root.predict(self,x)
    
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
        
    def __repr__(self):
        return self.__str__()

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

# 
# ## Expérimentations sur jeu de données artificielles 
# 
# + Prenez en main le code suivant, que génère-t-il comme données ? Visualiser quelques exemples.
# 
# + Sur quelques exemples, apprenez un arbre de décision et observer l'erreur. Tracer les frontières de décisions. 
# 
# + Observez l'erreur sur l'ensemble d'apprentissage. Comment se comporte-t-elle en fonction des deux paramètres ? Est-elle une bonne prédiction de l'erreur de votre modèle ? Comment de manière simple obtenir une meilleure prédiction (ensemble de test) ?
# 
# + Partager votre ensemble en deux sous-ensembles, un d'apprentissage qui vous servira à apprendre votre modèle, l'autre de test qui vous servira à évaluer l'erreur. Faites une évaluation intensive et tracez en fonction de la profondeur l'erreur en apprentissage et en test
# 
# + Ajouter un peu de bruit au données (paramètre epsilon de l'algorithme, compléter l'algo de manière a prendre en compte epsilon), recommencez vos expériences. Que remarquez-vous ?
# 

# <codecell>


def gen_arti(centerx=1,centery=1,sigma=0.1,nbex=1000,data_type=0,epsilon=0.02):
    #center : entre des gaussiennes
    #sigma : ecart type des gaussiennes
    #nbex : nombre d'exemples
    # ex_type : vrai pour gaussiennes, faux pour echiquier
    #epsilon : bruit

    if data_type==0:
        #melange de 2 gaussiennes
        xpos=np.random.multivariate_normal([centerx,centerx],np.diag([sigma,sigma]),nbex/2)
        xneg=np.random.multivariate_normal([-centerx,-centerx],np.diag([sigma,sigma]),nbex/2)
        data=np.vstack((xpos,xneg))
        y=np.hstack((np.ones(nbex/2),-np.ones(nbex/2)))
    if data_type==1:
        #melange de 4 gaussiennes
        xpos=np.vstack((np.random.multivariate_normal([centerx,centerx],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([-centerx,-centerx],np.diag([sigma,sigma]),nbex/4)))
        xneg=np.vstack((np.random.multivariate_normal([-centerx,centerx],np.diag([sigma,sigma]),nbex/4),np.random.multivariate_normal([centerx,-centerx],np.diag([sigma,sigma]),nbex/4)))
        data=np.vstack((xpos,xneg))
        y=np.hstack((np.ones(nbex/2),-np.ones(nbex/2)))

    if data_type==2:
        #echiquier
        data=np.reshape(np.random.uniform(-4,4,2*nbex),(nbex,2))
        y=np.ceil(data[:,0])+np.ceil(data[:,1])
        y=2*(y % 2)-1

    # un peu de bruit
    data[:,0]+=np.random.normal(0,epsilon,nbex)
    data[:,1]+=np.random.normal(0,epsilon,nbex)
    # on mélange les données
    idx = np.random.permutation((range(y.size)))
    data=data[idx,:]
    y=y[idx]
    return data,y

#affichage en 2D des donnees
def plot_data(x,labels):
    plt.scatter(x[labels<0,0],x[labels<0,1],c='red',marker='x')
    plt.scatter(x[labels>0,0],x[labels>0,1],c='green',marker='+')

#Frontiere de decision
def plot_frontiere(x,f,step=20): # script qui engendre une grille sur l'espace des exemples, calcule pour chaque point le label
                                # et trace la frontiere
    mmax=x.max(0)
    mmin=x.min(0)
    x1grid,x2grid=np.meshgrid(np.linspace(mmin[0],mmax[0],step),np.linspace(mmin[1],mmax[1],step))
    grid=np.hstack((x1grid.reshape(x1grid.size,1),x2grid.reshape(x2grid.size,1)))
    # calcul de la prediction pour chaque point de la grille
    res=np.array([f(grid[i,:]) for i in range(x1grid.size)])
    res=res.reshape(x1grid.shape)
    # tracer des frontieres
    plt.contourf(x1grid,x2grid,res,colors=('gray','blue'),levels=[-1,0,1])

# <markdowncell>

# ## Expérimentations sur USPS
# 
# Tester l'algorithme sur les données du [TME3 de MAPSI](http://webia.lip6.fr/~mapsi/pmwiki.php?n=Cours.TME3). Servez vous soit du code du TME3, soit du suivant pour lire le fichier de données. Attention ! L'implémentation est lourde, tester l'algorithme sur un arbre de petite profondeur. 
# 
# Visualisez sur quelques exemples les arbres. Choisissez dans la suite quelques classes parmi celles disponibles, partagez vos exemples en deux ensembles, un d'apprentissage et l'autre de test. Tracez l'erreur en apprentissage et en test en fonction de la profondeur et du nombre d'exemples en apprentissage. Que remarquez-vous ?

# <codecell>

def load_usps(filename):
    with open(filename,"r") as f:
        f.readline()
        data =[ [float(x) for x in l.split()] for l in f if len(l.split())>2]
    tmp = np.array(data)
    return tmp[:,1:],tmp[:,0].astype(int)

a,b = load_usps("2014_tme3_usps_train.txt")
print a[0],b[0]
exit()
# <markdowncell>

# ## Classification sur la base movielens 
# 
# ### Introduction
# 
# La base movielens est une base de données issue d'imdb, qui contient des informations sur des films (le genre, l'année de production, des tags) et des notes attribuées par les utilisateurs. Elle est utilisée généralement pour la recommendation de films. Nous allons l'utiliser dans le cadre de la classification, afin de prédire si un film est bon ou mauvais, dans deux contextes :
# 
# + en prenant en compte uniquement l'information sur le film et le score moyen du film
# 
# + en prenant en compte l'information de l'utilisateur qui score le film
# 
# Télécharger l'[archive suivante](http://www-connex.lip6.fr/~baskiotisn/Telecom/mvlenslight.zip)
# 
# Le bloc de code suivant est utilisé pour  charger et prétraiter les données.

# <codecell>

### TP arbres de decisions 
### Application aux donnees movielens : issues de imdb, donnees sur des films, des utilisateurs ayant notes les films
###
### Fichiers necessaires : 
### - users.dat : contient les informations sur les utilisateurs : 
###     * idUser, Sexe, Age, Profession, ZIP code
### - movies.dat : informations sur les films : 
###     * idMovie, Titre, [Genre]
### - ratings.dat : score des films :
###     *  idUser, idMovie, note (compris entre 1 et 5), timestamp
### - tags.dat : liste des tags
###     *  idTag, label, nombre de fois utilisé
### - tag_relevance.dat : adequation d'un tag a un film, entre 0 et 1
###     * idMovie, idTag, score
###


### modules :
### - numpy pour les outils mathématiques
### - matplotlib pour les outils graphiques

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import random

### liste des professions
occupation=[ "other", "academic/educator", "artist", "clerical/admin", "college/grad student", "customer service", "doctor/health care",
"executive/managerial", "farmer", "homemaker", "K-12 student", "lawyer", "programmer", "retired", "sales/marketing", "scientist", 
"self-employed", "technician/engineer", "tradesman/craftsman", "unemployed", "writer"]

### liste des genres
genre=['unknown','Action','Adventure','Animation',"Children's",'Comedy','Crime','Documentary','Drama','Fantasy','Film-Noir' ,
'Horror','Musical','Mystery','Romance','Sci-Fi','Thriller','War','Western']


### fonctions de lecture des fichiers 
def read_mlens(fname,sep='::'):
    """ Read a generic .dat file 
            - fname : file name
            - sep : separator 
        return : list of lists, each list contains the fields of each line read """
        
    def toint(x):
        if x.isdigit():
            return int(x)
        return x
    f=open(fname,'r')        
    tmp= [ s.lstrip().rstrip().split(sep) for s in f.readlines()]
    #lstrip et rstrip enleve les blancs de debut et de fin
    f.close()
    return [ [toint(x) for x in y] for y in tmp ]

def read_movies(fname="movies.dat"):
    """ Read movies information, reindexing movies
        return : res: binary matrix nbmovies x genre, res[i,j]= 1 if movie i has the genre j, 0 otherwise
                 movies2id : original idMovie to reindexed id
                 dic2Movies : reindexed id to list of movie information (id, title, genre)
    """                 
    movies=read_mlens(fname)
    dicMovies=dict()
    movies2id=dict()
    res=np.zeros((len(movies),len(genre)))
    for i,m in enumerate(movies):
        dicMovies[i]=m
        movies2id[m[0]]=i
        for g in m[2].split('|'):
            res[i,dicGenre[g]]=1
    return res,movies2id,dicMovies

def read_users(fname="users.dat"):
    """ Read users informations
        return : nbusers * 3 : gender (1 for M, 2 for F), age, occupation index"""
    users=read_mlens(fname)
    res=np.zeros((len(users),3))
    for u in users:
        res[u[0]-1,:]=[u[1]=='M' and 1 or 2, u[2],int(u[3])]
    return res
    
def read_files(mname="movies.dat",uname="users.dat",rname="ratings.dat"):
    """ Read all files 
        return :
            * movies: binary matrix movies x genre
            * users : matrix users x (gender, age, occupation index)
            * ratings : matrix movies x users, with score 1 to 5
            * movies2id : dictionary original id to reindexed id
            * dicMovies : dictionary reindexed id to movie information
    """
    print "Reading movies..."
    movies,movies2id,dicMovies=read_movies(mname)
    print "Reading users..."
    users=read_users(uname)
    print "Reading ratings..."
    rtmp=read_mlens(rname)
    ratings=np.zeros((movies.shape[0],users.shape[0]))
    for l in rtmp:
        ratings[movies2id[l[1]],l[0]-1]=l[2]        
    return movies,users,ratings,[],[],movies2id,dicMovies


def score2binary(score,thres=3):
    """transform score to labels 1/-1/0 : 1 if score >thres, -1 otherwise, 0 if no ratings"""
    return((score!=0)*(2*(score>thres)-1))

def get_movie_title(l,dicMovies):
    """ return : le titre des films contenues dans l """
    if (type(l)==type(1)):
        return [dicMovies[l][1]]
    return [dicMovies[x][1] for x in l]

def get_movie_year(l,dicMovies):
    """ return : liste des annees des films"""    
    if (type(l)==type(1)):
        return [dicMovies[l][1].split("(")[-1][:-1]]
    return [int(dicMovies[x][1].split("(")[-1][:-1]) for x in l]


def join_simple(users,movies,ratings):
    
    """ res : une matrice d'exemple, chaque ligne correspond a un utilisateur et un film qu'il a score :
        19 premieres colonnes : genre, 20 eme : annee, 3 dernieres : Sexe, Age, profession 
        ylab : 1 si score >3, 0 sinon
        yscore : le score de chaque rating
    """
    
    res=np.zeros((np.nonzero(ratings)[0].shape[0],users.shape[1]+1+movies.shape[1]))
    ylab=np.zeros(np.nonzero(ratings)[0].shape[0])
    yscore=np.zeros(np.nonzero(ratings)[0].shape[0])
    myear=get_movie_year(range(len(dicMovies)),dicMovies)
    cpt=0
    for m in range(ratings.shape[0]):
        for j in np.nonzero(ratings[m,:])[0]:  
            res[cpt,:]=np.hstack((movies[m,:],myear[m],users[j,:]))
            ylab[cpt]=score2binary(ratings[m,j])
            yscore[cpt]=ratings[m,j]
            cpt+=1
    return res,ylab,yscore

# <markdowncell>

# ### Prétraitement des données
# 
# Nous construisons d'abord la liste des genres des films (un film peut avoir plusieurs genres) et des professions des utilisateurs. Nous chargeons ensuite les données.

# <codecell>


### Preparation des donnees
### construction de la liste des genres et occupation
dicGenre=dict(zip(genre,range(len(genre))))
dicOccup=dict(zip(occupation,range(len(occupation))))

###lire les fichiers
movies,users,ratings,tagrelevance,tags,movies2id,dicMovies=read_files()

# <markdowncell>

# Les informations suivantes sont stockées :
# 
# + movies: une matrice binaire, chaque ligne un film, chaque colonne un genre, 1 indique le genre s'applique au film
# 
# + users : une matrice, chaque ligne un utilisateur, et les colonnes suivantes : sexe (1 masculin, 2 feminin),  age, index de la profession
# 
# + ratings : une matrice de score, chaque ligne un film, chaque colonne un utilisateur
# 
# + movies2id : dictionnaire permettant de faire la correspondance entre l'identifiant du film à l'identifiant réindexé
# 
# + dicMovies : dictionnaire inverse du précédent
# 
# ### Classification à partir de l'information unique du film
# 
# Notre matrice **movies** ne contient que les informations du genre du film. Nous voulons dans un premier temps y ajouter également l'année de production du film. Nous devons également constuire le score moyen du film. Les lignes suivantes le permettent.
# 
# + *<font style="BACKGROUND-COLOR: lightgray" color='red'> Sur quelques paramètres, que remarquez vous sur l'erreur d'apprentissage et de test ?</font>*
# 
# + *<font style="BACKGROUND-COLOR: lightgray" color='red'> La taille de l'ensemble de test joue-t-elle un rôle ?</font>*
# 
# + *<font style="BACKGROUND-COLOR: lightgray" color='red'> Tracer les courbes de ces deux erreurs en fonction de la profondeur. Que remarquez vous ? Quels sont les meilleurs paramètres pour l'erreur en apprentissage et en test ?</font>*
# 
# + *<font style="BACKGROUND-COLOR: lightgray" color='red'> Quelles sont les variables les plus importantes ?  </font>*
# 

# <codecell>

###moyenne des avis par film 
movieScore=np.nan_to_num(ratings.mean(1)*ratings.shape[1]/(1.*(ratings!=0).sum(1)))
movieScore[movieScore<=3]=-1
movieScore[movieScore>0]=1  

dbmovies=np.hstack((movies,np.reshape(get_movie_year(range(len(dicMovies)),dicMovies),(len(dicMovies),1))))

# <markdowncell>

# 
# ### Classification avec les informations utilisateurs
# 
# La fonction *join_simple* permet de fabriquer une matrice jointe où chaque ligne représente un film et un utilisateur. Les lignes de codes suivantes permettent de fabriquer la base d'apprentissage (*dataAll*) et les labels associés (*dataAllY*). Il est possible de passer une sous-matrice en paramètre (par exemple que les utilisateurs de moins de 20) et de fabriquer la base correspondante.
# 
# + *<font style="BACKGROUND-COLOR: lightgray" color='red'> Etudier comme précédement les erreurs en test et en apprentissage. Comparer ces erreurs si vous réduisez les utilisateurs par tranche d'age. Remarquez-vous des tranches d'age plus stable ?   </font>*

# <codecell>

dataAll,dataAllY,dataAllScore=join_simple(users,movies,ratings)

##utilisateurs de moins de 20 ans
uidx=(users[:,1]<20)
data20,data20Y,data20JScore=join_simple(users[uidx,:],movies,ratings[:,uidx])   

# <codecell>

