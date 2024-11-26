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
hidden_dim = 128
num_layers = 3
Animal = 0

model_params_filepath = 'trained_model_parameters_20241112.pth'

# feats_filepath = r'C:\Users\munib\Documents\Economo-Lab\data\ImputedMedullaDynamics\JPV8_2023-04-06_kinematics.csv'
feats_filepath = r'C:\Users\munib\Documents\Economo-Lab\data\ImputedMedullaDynamics\JPV11_2023-06-22_kinematics.csv'
trial_len = 515
num_trials = 213

# Load and preprocess features
feats = pd.read_csv(feats_filepath,header=None)
print(trial_len*num_trials, feats.shape[0])

X = preprocess_features(feats, trial_len, num_trials)

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


PCs = 8
Predicted_PCs, vExp = PCA(output, PCs, trial_len) # Returns Neural PCs and Plots Scree Plot

# After calculating Predicted_PCs, save it to a CSV file
np.savetxt("Predicted_PCs.csv", Predicted_PCs.reshape(trial_len*num_trials, PCs), delimiter=",")


# %%

fig, axs = plt.subplots(1, 2, figsize=(12, 8))

img1 = axs[0].imshow(Predicted_PCs[:, :, 0], aspect='auto', cmap='jet', origin='upper')
axs[0].set_xlabel('Time Points', fontweight='bold', fontsize=8)
axs[0].axvline(x=280, color='k', linestyle='--', linewidth=2)
axs[0].set_title('Predicted PC1')
plt.colorbar(img1, ax=axs[0], orientation='vertical')  

img2 = axs[1].imshow(Predicted_PCs[:, :, 1], aspect='auto', cmap='jet', origin='upper')
axs[1].set_xlabel('Time Points', fontweight='bold', fontsize=8)
axs[1].axvline(x=280, color='k', linestyle='--', linewidth=2)
axs[1].set_title('Predicted PC2')
plt.colorbar(img2, ax=axs[1], orientation='vertical') 

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.show()