#!env python

import pandas as pd
import numpy as np
import csv as csv
from sklearn.svm import SVC

# Data cleanup
# TRAIN DATA
train_df = pd.read_csv('cleanedTrain.csv', header=0)        # Load the train file into a dataframe
test_df  = pd.read_csv('cleanedTest.csv',  header=0)        # Load the test file into a dataframe

train_labels = train_df.Survived.values
test_id = test_df.PassengerId.values

train_df = train_df.drop(['Survived','PassengerId'], axis=1)
test_df = test_df.drop(['PassengerId'], axis=1)

# The data is now ready to go. So lets fit to the train, then predict to the test!
# Convert back to a numpy array
train_data = train_df.values
test_data = test_df.values

print 'Training...'
clf = SVC()
clf = clf.fit( train_data, train_labels )

print 'Predicting...'
output = clf.predict(test_data).astype(int)

predictions_file = open("CLF.csv", "wb")
open_file_object = csv.writer(predictions_file)
open_file_object.writerow(["PassengerId","Survived"])
open_file_object.writerows(zip(test_id, output))
predictions_file.close()
print 'Done.'

