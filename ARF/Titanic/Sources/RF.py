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


# female = 0, Male = 1
train_df['Gender'] = train_df['Sex'].map( {'female': 0, 'male': 1} ).astype(int)
test_df['Gender']  = test_df['Sex'].map( {'female': 0, 'male': 1} ).astype(int)

# Embarked from 'C', 'Q', 'S'
# All missing Embarked -> just make them embark from most common place
if len(train_df.Embarked[ train_df.Embarked.isnull() ]) > 0:
    train_df.Embarked[ train_df.Embarked.isnull() ] = train_df.Embarked.dropna().mode().values
if len(test_df.Embarked[ test_df.Embarked.isnull() ]) > 0:
    test_df.Embarked[ test_df.Embarked.isnull() ] = test_df.Embarked.dropna().mode().values
# convert all Embarked strings to int
Ports = list(enumerate(np.unique(train_df['Embarked']))) # determine all values of Embarked,
Ports_dict = { name : i for i, name in Ports }           # set up a dictionary in the form  Ports : index
train_df.Embarked = train_df.Embarked.map( lambda x: Ports_dict[x]).astype(int) # Convert all Embark strings to int
test_df.Embarked = test_df.Embarked.map( lambda x: Ports_dict[x]).astype(int)


# Remove the Name column, Cabin, Ticket
train_df = train_df.drop(['Survived','Name', 'Ticket', 'Cabin', 'PassengerId'], axis=1)
test_df  = test_df.drop(['Survived','Name', 'Ticket', 'Cabin', 'PassengerId'], axis=1)

# The data is now ready to go. So lets fit to the train, then predict to the test!
# Convert back to a numpy array
train_data = train_df.values
test_data = test_df.values

print 'Training...'
forest = RandomForestClassifier(n_estimators=100)
forest = forest.fit( train_data, train_labels )

print 'Predicting...'
output = forest.predict(test_data).astype(int)
print output
exit()
predictions_file = open("myfirstforest.csv", "wb")
open_file_object = csv.writer(predictions_file)
open_file_object.writerow(["PassengerId","Survived"])
open_file_object.writerows(zip(ids, output))
predictions_file.close()
print 'Done.'
