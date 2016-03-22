# Import the random forest package
from sklearn.ensemble import RandomForestClassifier 
from sklearn import cross_validation
import numpy as np
dataset = np.loadtxt('training_data.csv', delimiter=",")

# Create the random forest object which will include all the parameters
# for the fit
forest = RandomForestClassifier(n_estimators = 100)

# Fit the training data to the Survived labels and create the decision trees
forest_fit = forest.fit(dataset[0::,1::],dataset[0::,0])

importances = forest.feature_importances_
std = np.std([tree.feature_importances_ for tree in forest.estimators_],
             axis=0)

# Take the same decision trees and run it on the test data
#output = forest_fit.predict(test_data)

#cross cross_validation
scores = cross_validation.cross_val_score(forest_fit, dataset[0::,1::], dataset[0::,0], cv=5)

s = pickle.dump(forest_fit)