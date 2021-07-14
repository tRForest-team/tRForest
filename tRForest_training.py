#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Import relevant libraries
import numpy as np
from pandas import read_csv
import csv
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.metrics import confusion_matrix
from sklearn.metrics import classification_report
from scipy import stats
import matplotlib.pyplot as plt
from sklearn.metrics import plot_roc_curve
from sklearn.model_selection import RandomizedSearchCV
from random import randint
import joblib

# Define random forest classifier function
def RFC(arr, num_trees, max_features, ftrs):
    X = arr[:,ftrs]
    Y = arr[:,0] 
    kfold = KFold(n_splits=10, random_state=7)
    model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features)
    results = cross_val_score(model, X, Y, cv=kfold)
    return results.mean()

# Get feature data and convert to an array of values
filename = '/path/to/file/rf_features.csv'
names = ['trf_id', 'trfdb_id', 'trf_sequence', 'fixed_trf_sequence', 'mrna_id', 'enst_id', 'mrna_transcript', 'fixed_mrna_transcript', 'chr', 'start_loc', 'end_loc', 'length_utr', 'trf_binding_loc', 'mrna_binding_loc', 'binding_energy', 'seed', 'au_content', 'num_paired_pos', 'binding_region_length', 'longest_consecutive', 'pos_longest_consecutive', 'three_prime_pairs', 'seed_end_diff', 'accessibility', 'me_motif', 'phylop_stem', 'phylop_flanking', 'phastcons_stem', 'phastcons_flanking', 'istarget']
dataframe = read_csv(filename, names=names, low_memory=False)
array = dataframe.values

# Set hyperparameters and split up array into features and true value for training

num_trees = 200
max_features = 'auto'
max_depth = None
min_samples_split = 2
min_samples_leaf = 1
bootstrap = False

X = array[:, [14,15,16,17,18,19,20,21,22,25,26]]
X = X.astype(float)
Y = array[:,29]
Y = Y.astype(float)

# Training with 10-fold cross-validation
kfold = KFold(n_splits=10, random_state=7, shuffle=True)
model = RandomForestClassifier(n_estimators=num_trees, max_features=max_features, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, bootstrap=bootstrap) 
scoring = 'roc_auc'
results = cross_val_score(model, X, Y, cv=kfold, scoring=scoring)
accuracy_res = cross_val_score(model, X, Y, cv=kfold, scoring='accuracy')

# Training with classic 67%:33% split of training and testing data
test_size = 0.33
seed = 7
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=seed)
model_split = RandomForestClassifier(n_estimators=num_trees, max_features=max_features, max_depth=max_depth, min_samples_split=min_samples_split, min_samples_leaf=min_samples_leaf, bootstrap=bootstrap)
model_split.fit(X_train, Y_train) 
predicted = model_split.predict(X_test)

# Fit model

model.fit(X_train, Y_train)

# Obtain metrics
matrix = confusion_matrix(Y_test, predicted)
report = classification_report(Y_test, predicted, output_dict=True)
ax = plt.gca()
rfc_disp = plot_roc_curve(model_split, X_test, Y_test, ax=ax, alpha=0.8)
plt.title("Receiver Operating Characteristic Curve")
plt.show()

print(" Accuracy: " + str(round(accuracy_res.mean(), 4)) + " (" + str(round(stats.sem(accuracy_res), 4)) + ")")
print("      AUC: " + str(round(results.mean(), 4)) + " (" + str(round(stats.sem(results), 4)) + ")")

precision_0 = list((list(report.items())[0][1]).items())[0][1]
recall_0 = list((list(report.items())[0][1]).items())[1][1]
precision_1 = list((list(report.items())[1][1]).items())[0][1]
recall_1 = list((list(report.items())[1][1]).items())[1][1]
f1score_0 = list((list(report.items())[0][1]).items())[2][1]
f1score_1 = list((list(report.items())[1][1]).items())[2][1]
print("Precision: " + str(round(np.mean(np.array([precision_0, precision_1])), 4)))
print("   Recall: " + str(round(np.mean(np.array([recall_0, recall_1])), 4)))
print(" F1 Score: " + str(round(np.mean(np.array([f1score_0, f1score_1])), 4)))

# Save models
joblib.dump(model, "/Users/RohanParikh/Google Drive/Dutta Lab/random forests/tRForest_kfold.joblib")
joblib.dump(model_split, "/Users/RohanParikh/Google Drive/Dutta Lab/random forests/tRForest_split.joblib")

