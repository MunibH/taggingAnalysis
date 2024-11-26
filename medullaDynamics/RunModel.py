# %%
# conda activate ephys-torch

# %%

import sys
import pandas as pd
import torch
from model import LSTM_model
from preprocessing import preprocess_features
import matplotlib.pyplot as plt 
import numpy as np

# %%

# Command-line arguments
savepth = sys.argv[1]
hidden_dim = int(sys.argv[2])
num_layers = int(sys.argv[3])
Animal = int(sys.argv[4])
trial_len = int(sys.argv[5])
model_params_filepath = sys.argv[6]
feats_filepath = sys.argv[7]
PCs = int(sys.argv[8])
num_trials = int(sys.argv[9])

# Load and preprocess features
feats = pd.read_csv(feats_filepath,header=None) # loads kinematc data of size (nTime*nTrials,nFeats)
X = preprocess_features(feats, trial_len, num_trials) # returns nan-imputed, zscored kinematic data of size (nTrials,nTime,nFeats)

# Initialize and load model
input_dim = X.shape[2]
device = torch.device('cpu')
model = LSTM_model(input_dim, hidden_dim, num_layers, Animal)
model.load_state_dict(torch.load(model_params_filepath, map_location=device))

# Get model output
with torch.no_grad():
    output = model(X)


# %%

# PCA on LSTM OUTPUT 
def PCA(output, PCs, trial_len): 
    Y = output.cpu().numpy()

    # Reshape to (time points * trials, 128 Dims)
    Y = Y.reshape((Y.shape[0] * Y.shape[1], Y.shape[2]))
    
    # Perform SVD
    A = (Y - np.mean(Y, axis=0)).T
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt.T

    # Cumulative Variance 
    vExp = [S[i]/np.sum(S) for i in range(len(S))]
    cVar = []

    temp_var = 0.0
    for i in range(len(S)):
        temp_var += vExp[i]
        cVar.append(temp_var)

    # # Scree Plot
    # plt.plot(range(1, 1 + len(S)),cVar, color = 'k')
    # plt.ylim((0,1))
    # plt.xlim((1, len(S)))
    # plt.xlabel('PC')
    # plt.ylabel('Cumulative Variance')
    # plt.axvline(x = PCs, color = 'r', linestyle = '--', label = f'PC{PCs} - vExp = {cVar[PCs]:.3f}')
    # plt.title('Scree Plot of LSTM Dimensions', fontweight = 'bold')
    # plt.legend()

    # Perform Dimensionality Reduction
    S = np.diag(S)
    Y_pc = V[:, :PCs] @ S[:PCs, :PCs]
    Y_pc = Y_pc.reshape((Y_pc.shape[0] // trial_len, trial_len, PCs))


    return Y_pc, vExp[:PCs]

Predicted_PCs, vExp = PCA(output, PCs, trial_len) # Returns Neural PCs and Plots Scree Plot
# Predicted_Pcs = Predicted_PCs.transpose(Predicted_PCs, (1, 0, 2))

    # # Reshape X from (time*trials, feats) to (trials, time, feats)
    # X_reshaped = X.reshape((trial_len, num_trials, X.shape[1])).astype(np.float32)
    # X_reshaped = np.transpose(X_reshaped, (1, 0, 2))

# After calculating Predicted_PCs, save it to a CSV file
np.savetxt(savepth + "_PCs.csv", Predicted_PCs.reshape(trial_len*num_trials, PCs), delimiter=",")
np.savetxt(savepth + "_vExp.csv", vExp, delimiter=",")
