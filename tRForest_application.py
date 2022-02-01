#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Load models
model = joblib.load("/path/to/tRForest_kfold.joblib")
model_split = joblib.load("/path/to/tRForest_split.joblib")

# Set paths and get filenames in relevant directory
from os import walk
mypath = "/path/to/all_feature_profiles/" # Directory with feature profiles of all tRFs
outputpath = "/path/to/output/" # Output directory

f = []
for (dirpath, dirnames, filenames) in walk(mypath):
    f.extend(filenames)
    break
    
# Generate predictions for feature profiles
for k in range(len(f)):
    new_dataset_filename = mypath + f[k]
    names = ['trfdb_id','trf_sequence','tran_id','tran_ver','name','mrna_sequence','chr','start_loc','end_loc','length_utr','trf_binding_loc', 'mrna_binding_loc', 'binding_energy', 'seed', 'au_content', 'num_paired_pos', 'binding_region_length', 'longest_consecutive', 'pos_longest_consecutive', 'three_prime_pairs', 'seed_end_diff', 'phylop_stem', 'phylop_flanking']
    
    # Get feature values
    X_new = read_csv(new_dataset_filename, names=names, low_memory=False)
    if(len(X_new) == 2):
        X_new = X_new[1:2]
    else:
        X_new = X_new[1:(len(X_new)-1)]
    X_new = X_new.values
    X_new_features = X_new[:, 12:23]
    X_new_features = X_new_features.astype(float)
    X_new_features = np.nan_to_num(X_new_features)
    
    # Get predictions with and without probabilities
    Y_prednew_kfold = model.predict(X_new_features)
    Y_predproba_kfold = model.predict_proba(X_new_features)[:,1]
    
    # Output files
    idx_start = f[k].index('files_')
    idx_end = f[k].index('.csv')

    resultFile = outputpath + 'kfold_preds/trf-' + f[k][idx_start+6:idx_end] + '_kfold_pred.csv'
    wtr = csv.writer(open (resultFile, 'w'), delimiter=' ', lineterminator='\n')
    for x in Y_prednew_kfold : wtr.writerow ([x])
    
    resultFile2 = outputpath + 'kfold_preds_proba/trf-' + f[k][idx_start+6:idx_end] + '_kfold_predproba.csv'
    wtr = csv.writer(open (resultFile2, 'w'), delimiter=' ', lineterminator='\n')
    for x in Y_predproba_kfold : wtr.writerow ([x])

