# %%
# conda activate facemap

import os
import torch
from train import TrainFacemapNetwork

# %%

datapth = r'C:\Users\munib\Documents\Economo-Lab\data'
proj = 'Facemap'
datapth = os.path.join(datapth,proj)

metadata = open('metadata.csv')
files = [i.rstrip().replace(" ", "_") + "_FacemapData.mat" for i in metadata]

# %% MODEL PARAMETERS

modelpar = {
    'delay': -1,  # negative, neural activity precedes behavior, positive, succeeds
    'smoothing_penalty': 0,  # weighting factor for weight change regularization
    'n_iter': 150,  # num epochs
    'learning_rate': 1e-3,
    'annealing_steps': 1,  # number of annealing steps of the learning rate
    'weight_decay': 1e-4,  # L2 regularizer
    'n_filt': 10,
    'kernel_size': 201,
    'n_core_layers': 2,
    'n_latents': 500,
    'n_out_layers': 1,
    'n_med': 50,
    'relu_wavelets': True,  # True is good
    'relu_latents': True  # True is good
}

device = torch.device("cuda")
# device = torch.device("cpu")

# %%

for f in files:
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('---------Training facemap model - ' + f.split('_FacemapData.m')[0] + '---------')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    TrainFacemapNetwork(datapth,f,modelpar,device)
    
    
    
    

# %%
