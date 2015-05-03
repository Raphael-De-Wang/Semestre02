#!env python

import pandas as pd
import numpy as np
import csv as csv
from sklearn.ensemble import RandomForestClassifier

# Data cleanup
# TRAIN DATA
train_df = pd.read_csv('cleanedTrain.csv', header=0)        # Load the train file into a dataframe
test_df  = pd.read_csv('cleanedTest.csv',  header=0)        # Load the test file into a dataframe

train_labels = train_df.Survived.values
test_id = test_df.PassengerId.values

train_df = train_df[['AgeClass', 'Age', 'Fare', 'Fare_Per_Person', 'Family_Size', 'Embarked', 'Deck', 'SibSp', 'Title','Gender']]
test_df = test_df[['AgeClass', 'Age', 'Fare', 'Fare_Per_Person', 'Family_Size', 'Embarked', 'Deck', 'SibSp', 'Title','Gender']]

# The data is now ready to go. So lets fit to the train, then predict to the test!
# Convert back to a numpy array
train_data = train_df.values
test_data = test_df.values

print 'Training...'
forest = RandomForestClassifier(n_estimators=500,max_depth=3)
forest = forest.fit( train_data, train_labels )

print 'Predicting...'
output = forest.predict(test_data).astype(int)

predictions_file = open("RF.csv", "wb")
open_file_object = csv.writer(predictions_file)
open_file_object.writerow(["PassengerId","Survived"])
open_file_object.writerows(zip(test_id, output))
predictions_file.close()
print 'Done.'

