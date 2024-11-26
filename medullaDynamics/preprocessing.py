# Preprocessing

# Import Libraries 
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
import torch 
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import zscore
import random

def preprocess_features(X, trial_len, num_trials):
    # Drop 'Trial Type' column if it exists
    if 'Trial Type' in X.columns:
        X = X.drop(columns=['Trial Type'])
    
    # Convert to numpy array
    X = X.to_numpy()
    
    # Impute missing values with the column mean
    col_means = np.nanmean(X, axis=0)
    inds = np.where(np.isnan(X))
    X[inds] = np.take(col_means, inds[1])
    
    # Z-score the kinematic features
    X = zscore(X, axis=0)

    # # Truncate X to fit an integer number of trials (if necessary)
    # X = X[:n_trials * trial_len, :]

    # Reshape X from (time*trials, feats) to (trials, time, feats)
    # X_reshaped = X.reshape((trial_len, num_trials, X.shape[1])).astype(np.float32)
    # X_reshaped = np.transpose(X_reshaped, (1, 0, 2))
    X_reshaped = X.reshape((num_trials, trial_len, X.shape[1])).astype(np.float32)
    X_reshaped = torch.from_numpy(X_reshaped)
    
    return X_reshaped


# def preprocess_features(X, trial_len):

#     # Drop Trial Type Column if exists 
#     if 'Trial Type' in X.columns: 
#         X = X.drop(columns = ['Trial Type'])
    
#     X = X.to_numpy()
    
#     # Z - scoring Kinematic Features (X)
#     X = zscore(X, axis = 0)

#     # Reshape (Trials, Time Points, Features + ID)
    
    
#     # X_reshaped = X.reshape((X.shape[0] // trial_len, trial_len, X.shape[1])).astype(np.float32)
#     # X_reshaped = torch.from_numpy(X_reshaped)
    
#     # Determine the number of trials
#     n_trials = X.shape[0] // trial_len

#     # Truncate X to fit an integer number of trials (if necessary)
#     X = X[:n_trials * trial_len, :]

#     # Reshape X from (time*trials, feats) to (trials, time, feats)
#     X_reshaped = X.reshape((n_trials, trial_len, X.shape[1])).astype(np.float32)
#     X_reshaped = torch.from_numpy(X_reshaped)


#     return X_reshaped