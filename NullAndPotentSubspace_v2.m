clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'manopt')));
addpath(genpath(fullfile(utilspth,'subspace')));

clc

%% TODO
% 1) dimensionality of each subspace...
% 2) inpar.delay and inpar.responseOnly should be more robust. I just hardcoded
%       the time values right now
% 3) deal with dual probes for visualization and analysis after finding
%       subspaces

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

params.behav_only = 0;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);
meta1 = ALM_SessionMeta(meta,datapth);
meta2 = tjM1_SessionMeta(meta,datapth);
meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth);  % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

meta = meta(1:2);

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'regress'; % 'st' or 'ta' or '2pca' or 'regress'

% inpar.trials = 'all'; % specify 'all' or condition numbers
inpar.trials = 'all'; % 2:5;

inpar.delayOnly = false; % if true, only use delay epoch for subspace estimation
inpar.responseOnly = false; % if true, only use response epoch for subspace estimation

% dimensionality (will soon change to dynamically set this)
inpar.nNullDim = 4;
inpar.nPotentDim = 4;

inpar.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)

inpar.estimateDimensionality = false; % if true, find upper-bound of dimensionality using parallel analysis
inpar.dimpth = 'results\Dimensionality'; % where to save parallel analysis results

inpar.standardize = 'baseline'; % 'baseline' or 'trial' - which time points to compute zscore stats over

%% LOAD DATA and PERFORM SUBSPACE ANALYSIS

all_alignment = struct();
all_ve.null = nan(size(meta));
all_ve.potent = nan(size(meta));
all_contrib = struct();
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth all_alignment all_ve all_contrib
    
    % load data
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta)) ' : ' meta(isess).anm ' ' meta(isess).date])
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);

    % 
    SubspaceEngine(thismeta,inpar,sessobj,sesspar,params,tag,me);

end