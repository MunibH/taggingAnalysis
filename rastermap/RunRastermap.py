# %%
# conda activate rastermap

# %%
from scipy.io import loadmat

import numpy as np
import matplotlib.pyplot as plt
# importing rastermap
# (this will be slow the first time since it is compiling the numba functions)
from rastermap import Rastermap, utils
from scipy.stats import zscore

import os

from sklearn.preprocessing import StandardScaler

# %%

filepath = r'C:\Users\munib\Documents\Economo-Lab\data\rastermap'
filename = 'MAH24_2024-06-11_RasterMapData.mat'
data = loadmat(os.path.join(filepath,filename))
# nTimeEachTrial = data['nTimeEachTrial']
# motionSVD = data['videoSVD']
# neuralActivity = data['neuralActivity']
# dt = data['dt'][0][0]

anm = str(data['anm'][0])
session_date = str(data['sessiondate'][0])
spks = data['spks'] # binned firing rates of size (neurons x time bins)
trial_start = data['trialstart']
sample = data['sample']
delay = data['delay']
gocue = data['gocue']
dt = data['dt'][0]
tedges = data['tedges']
tagged = data['tagged']
psth = data['psth']
time = data['time']

# %%

model = Rastermap(n_clusters=6, # None turns off clustering and sorts single neurons 
                  n_PCs=20, # use fewer PCs than neurons
                  locality=0.1, # some locality in sorting (this is a value from 0-1)
                  time_lag_window=15, # use future timepoints to compute correlation
                  grid_upsample=1, # 0 turns off upsampling since we're using single neurons
                  mean_time=False, # project out mean across neurons at each timepoint
                ).fit(spks)
y = model.embedding # neurons x 1
isort = model.isort

assigned_cluster = model.embedding_clust
cluster_activity = model.X_nodes
neuron_pc_embedding = model.Usv
cluster_pc_embedding = model.U_nodes

# %%

# # for each assigned cluster, plot mean cluster activity, and any tagged units separately
# cluster_ids = np.unique(assigned_cluster)
# for iclu in cluster_ids:
#     cluster_mask = assigned_cluster==iclu
#     cluster_mean = np.mean(psth[:,cluster_mask,:],axis=1)
    
#     f, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2), gridspec_kw={'width_ratios': [1, 1]})    
#     ax1.plot(time,cluster_mean[:,0],c='b',lw=1.5)
#     ax1.plot(time,cluster_mean[:,1],c='r',lw=1.5)
#     ax1.set_title(f'Cluster {iclu} Mean, {np.sum(cluster_mask)} units')
#     ax1.set_xlabel("Time from go cue (s)")
#     ax1.set_ylabel("spks/sec")
    
#     istag = np.where((tagged.flatten() == 1) & (assigned_cluster.flatten() == iclu))[0]
#     if istag.size!=0:
#         for itag in istag:
#             tag_mean = (psth[:,itag,:])
#             ax2.plot(time, tag_mean[:, 0], c='b', lw=1)
#             ax2.plot(time, tag_mean[:, 1], c='r', lw=1)
#             ax2.set_title(f'Tagged Units in Cluster {iclu}')
            
    
#     plt.show()
    
    
# for each assigned cluster, plot mean cluster activity, and any tagged units separately
cluster_ids = np.unique(assigned_cluster)
for iclu in cluster_ids:
    cluster_mask = assigned_cluster==iclu
    cluster_mean = np.mean(psth[:,cluster_mask,:],axis=(1,2))
    
    f, ax = plt.subplots(1, 1, figsize=(4, 3))    
    ax.plot(time,cluster_mean,c=(0.3, 0.3, 0.3),lw=1.8)

    
    istag = np.where((tagged.flatten() == 1) & (assigned_cluster.flatten() == iclu))[0]
    if istag.size!=0:
        for itag in istag:
            tag_mean = np.mean(psth[:,itag,:],axis=1)
            ax.plot(time, tag_mean, c=(0.9882, 0.4667, 0.0118), lw=1.3)
    yl = ax.get_ylim()
    ax.plot(np.array([0,0]),yl,c='k',ls='--')
    ax.plot(np.array([-1.2,-1.2]),yl,c='k',ls='--')
    ax.plot(np.array([-1.85,-1.85]),yl,c='k',ls='--')
    ax.set_ylim(yl)
    ax.set_title(f'Cluster{iclu}, {np.sum(cluster_mask)} units', fontsize=12, pad=1)
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.tick_params(axis='both', which='major', labelsize=12,length=7,width=2)
    ax.set_xlabel("Time from go cue (s)", fontsize=14)
    ax.set_ylabel("spks/sec", fontsize=14)
    plt.tight_layout()
    plt.show()

    

# %% cluster activity

# timepoints to visualize
xmin = 70000
xmax = xmin + 3000

# make figure with grid for easy plotting
fig = plt.figure(figsize=(12,8), dpi=200)
grid = plt.GridSpec(10, 24, figure=fig, wspace = 0.1, hspace = 0.4)

# plot sorted neural activity
ax = plt.subplot(grid[2:, :-5])
ax.imshow(cluster_activity[:, xmin:xmax], cmap="inferno",aspect='auto')
ax.set_xlabel("time")
ax.set_ylabel("clusters")

# %% PLOT NEURONS SORTED BY ISORT

# timepoints to visualize
xmin = 70000
xmax = xmin + 4000

# make figure with grid for easy plotting
fig = plt.figure(figsize=(12,8), dpi=200)
grid = plt.GridSpec(10, 24, figure=fig, wspace = 0.1, hspace = 0.4)


# plot sorted neural activity
ax = plt.subplot(grid[2:, :-5])
ax.imshow(spks[isort, xmin:xmax], cmap="inferno",vmin=0,vmax=100,aspect='auto')
# ax.imshow(spks[isort, xmin:xmax], cmap="inferno", vmin=0, vmax=1.2, aspect="auto")
ax.set_xlabel("time")
ax.set_ylabel("neurons")

# excitatory cells in yellow, and inhibitory cells in dark blue
# (could replace this with a colorbar or other property)
ax = plt.subplot(grid[2:, -5])
ax.imshow(tagged[isort, np.newaxis],
          cmap="viridis", aspect="auto")
ax.axis("off")

# %% CREATE SUPERNEURONS FROM RASTERMAP, SORT DATA ND SUM OVER NEIGHBORING NEURONS

nbin = 4 # number of neurons to bin over 
sn = utils.bin1d(spks[isort], bin_size=nbin, axis=0) # bin over neuron axis

# timepoints to visualize
xmin = 70000
xmax = xmin + 2000

# make figure with grid for easy plotting
fig = plt.figure(figsize=(12,8), dpi=200)
grid = plt.GridSpec(10, 24, figure=fig, wspace = 0.1, hspace = 0.4)


# plot sorted neural activity
ax = plt.subplot(grid[2:, :-5])
ax.imshow(sn[:, xmin:xmax], cmap="inferno",vmin=0,vmax=100,aspect='auto')
# ax.imshow(spks[isort, xmin:xmax], cmap="inferno", vmin=0, vmax=1.2, aspect="auto")
ax.set_xlabel("time")
ax.set_ylabel("superneurons")


# %%

