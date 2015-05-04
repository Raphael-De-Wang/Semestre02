#!env python

import pandas as pd
import numpy as np
import csv as csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.feature_selection import RFECV
from sklearn.ensemble import RandomForestClassifier


train_df = pd.read_csv('cleanedTrain.csv', header=0)        # Load the train file into a dataframe

features_list = train_df.columns.values[1::]
y = train_df[["PassengerId","Survived"]].values[:, 1]
train_df = train_df.drop(["PassengerId","Survived"],axis=1)
X = train_df.values

# Fit a random forest with (mostly) default parameters to determine feature importance
forest = RandomForestClassifier(oob_score=True, n_estimators=200)
selector = RFECV(forest, step=1, cv=5)
selector = selector.fit(X, y)
print features_list
print 'selector.support_ : ', selector.support_
print 'selector.ranking_ : ', selector.ranking_
print 'selector.grid_scores_ : ', selector.grid_scores_
print 'selector.n_features_  : ', selector.n_features_ 
