# %%

import torch
from model import KeypointsNetwork

from utils import LoadFacemapData

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import zscore

import os
import numpy as np

from scipy.io import savemat

# %%

def TrainFacemapNetwork(datapth,filename,modelpar,device):
    

    # %% parameters

    delay = modelpar['delay']
    smoothing_penalty = modelpar['smoothing_penalty']
    n_iter = modelpar['n_iter']
    learning_rate = modelpar['learning_rate']
    annealing_steps = modelpar['annealing_steps']
    weight_decay = modelpar['weight_decay']

    n_filt = modelpar['n_filt']
    kernel_size = modelpar['kernel_size']
    n_core_layers = modelpar['n_core_layers']
    n_latents = modelpar['n_latents']
    n_out_layers = modelpar['n_out_layers']
    n_med = modelpar['n_med']
    relu_wavelets = modelpar['relu_wavelets']
    relu_latents = modelpar['relu_latents']


    # %% get data, zscore

    nTimeEachTrial, videoSVD, neuralActivity, dt = LoadFacemapData(filename,datapth)
    # nTimeEachTrial - index of each trial end in first axis of motionSVD and neuralActivity
    #   size (1,trials)
    # videoSVD - size(time*trials,100 motSVs + 100 movSVs)
    # neuralActivity - size(time*trials, neurons)

    motionSVD = zscore(videoSVD)
    neuralActivity = neuralActivity - np.mean(neuralActivity,0)
    nTime,nSV = motionSVD.shape
    _,nUnits = neuralActivity.shape
    tcam = np.arange(nTime)
    tneural = np.arange(nTime)

    # %% instantiate and train model

    np.random.seed(0)
    torch.manual_seed(0)
    torch.cuda.manual_seed(0)

    # instantiate model
    model = KeypointsNetwork(
        n_in=motionSVD.shape[-1], n_out=neuralActivity.shape[-1],
        n_filt=n_filt,kernel_size=kernel_size,n_core_layers=n_core_layers,
        n_latents=n_latents,n_out_layers=n_out_layers,n_med=n_med,
        relu_wavelets=relu_wavelets,relu_latents=relu_latents,identity=False
    ).to(device)

    # train model
    (
        y_pred_test,
        varexps_identity,
        spks_pred_test,
        varexps_identity_neurons,
        itest, epoch_train_loss, epoch_lr, epoch_ve
    ) = model.train_model(motionSVD, neuralActivity, tcam, tneural, 
                        smoothing_penalty=smoothing_penalty, delay=delay,weight_decay=weight_decay,
                        n_iter=n_iter,learning_rate=learning_rate,annealing_steps=annealing_steps,
                        device=device)

    # Save model weights
    savepth = datapth
    savename = filename.split('_FacemapData.m')[0]
    torch.save(model.state_dict(), os.path.join(savepth,savename) + '_ModelWeights.pth')
    # model.load_state_dict(torch.load('deep_network_model.pth'))

    # %% test 

    X_test = torch.tensor(motionSVD.copy(), dtype=torch.float32).unsqueeze(0).to(device) # (batch,time,feats)

    with torch.no_grad():  # Disable gradient computation for inference
        predictions, deep_behavioral_features = model(X_test)
    # # forward(self, x, sample_inds=None, animal_id=0)
    # # torch.Size([1, 5741, 500]) X_train[animal][batch]

    predictions = predictions.cpu().detach().numpy().squeeze() # (time,units)
    
    # %% save
    mdic = {"pred" : predictions, 
        "neuralActivity" : neuralActivity,
        "motionSVD" : motionSVD,
        "nTimeEachTrial" : nTimeEachTrial,
        "dt" : dt,
        "ve_neuron" : varexps_identity_neurons,
        "epoch_train_loss" : epoch_train_loss,
        "epoch_lr" : epoch_lr,
        "epoch_ve" : epoch_ve}

    savepth = datapth
    savename = filename.split('_FacemapData.m')[0]
    savemat(os.path.join(savepth,savename) + '_FacemapPredictions.mat' , mdic)
