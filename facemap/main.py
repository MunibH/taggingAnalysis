# %%
# train.py
import torch
from model import KeypointsNetwork

from utils import LoadFacemapData

import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import zscore

# %% get data, zscore

fpth = r'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data'
fname = 'JPV11_2023-06-16_FacemapData.mat'
nTimeEachTrial, videoSVD, neuralActivity, dt = LoadFacemapData(fpth,fname)
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

delay = -1 # negative, neural activity precedes behavior, positive, succeeds
smoothing_penalty = 0 # weighting factor for weight change regularization
n_iter = 150 # num epochs
learning_rate = 1e-3
annealing_steps = 1 # number of annealing steps of the learning rate
weight_decay = 1e-4 # L2 regularizer

n_filt=10
kernel_size=201
n_core_layers=2
n_latents=500
n_out_layers=1
n_med=50
relu_wavelets=True # True is good
relu_latents=True # True is good

device = torch.device("cuda")
# device = torch.device("cpu")

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

# After training, you can save the model if needed
torch.save(model.state_dict(), 'deep_network_model.pth')
# model.load_state_dict(torch.load('deep_network_model.pth'))

# %% plot training

plotve = epoch_ve.copy()
plotve[plotve==0] = np.nan
indices_not_nan = np.where(~np.isnan(plotve))[0]

f = plt.figure()
ax = plt.gca()
plt.plot(np.arange(n_iter),epoch_train_loss)
plt.scatter(np.arange(n_iter),epoch_train_loss)
plt.plot(indices_not_nan,plotve[indices_not_nan])
plt.scatter(indices_not_nan,plotve[indices_not_nan])
plt.xlabel('Epoch')
plt.ylabel('Validation loss')
plt.show()

f = plt.figure()
ax = plt.gca()
plt.plot(np.arange(n_iter),epoch_lr)
plt.scatter(np.arange(n_iter),epoch_lr)
plt.xlabel('Epoch')
plt.ylabel('Learning rate')
ax.set_yscale('log')
plt.show()

# %% test 

X_test = torch.tensor(motionSVD.copy(), dtype=torch.float32).unsqueeze(0).to(device) # (batch,time,feats)

with torch.no_grad():  # Disable gradient computation for inference
    predictions, deep_behavioral_features = model(X_test)
# # forward(self, x, sample_inds=None, animal_id=0)
# # torch.Size([1, 5741, 500]) X_train[animal][batch]

predictions = predictions.cpu().detach().numpy().squeeze() # (time,units)



# %% plot predictions

nTimes = np.cumsum(nTimeEachTrial)

trial = 100
unit = 52

ii = np.where(varexps_identity_neurons<0)[0]
# plt.scatter(np.arange(len(varexps_identity_neurons)),varexps_identity_neurons)

if trial==0:
    itime = np.arange(0,nTimes[trial])
else:
    itime = np.arange(nTimes[trial-1],nTimes[trial])
    
    
f = plt.figure
ax = plt.gca
plt.plot(neuralActivity[itime,unit])
plt.plot(predictions[itime,unit])
plt.title(varexps_identity_neurons[unit])
plt.show()


# %% save predictions

from scipy.io import savemat
import os
mdic = {"pred" : predictions, 
        "neuralActivity" : neuralActivity,
        "motionSVD" : motionSVD,
        "nTimeEachTrial" : nTimeEachTrial,
        "dt" : dt,
        "ve_neuron" : varexps_identity_neurons,
        "epoch_train_loss" : epoch_train_loss,
        "epoch_lr" : epoch_lr,
        "epoch_ve" : epoch_ve}
# fname = 'MAH24_2024-06-11_FacemapPredictions.mat'
# fname = 'MAH24_2024-06-11_FacemapPredictions_noAnneal.mat'
fname = 'JPV11_2023-06-16_FacemapPredictions.mat'
savemat(os.path.join('data',fname) , mdic)
