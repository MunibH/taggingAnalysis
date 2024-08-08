clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warp params
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
% 2afc
params.condition(1)     = {'R&hit&~stim.enable&~autowater'};                 % right hits, no stim, aw off   (1)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off    (2)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off  (3)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off   (4)

% for projections
% 2afc
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off   (5)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off    (6)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off  (7)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off   (8)


% for ramping
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all hits, no stim, aw off (9)

% wc
params.condition(end+1) = {'R&hit&~stim.enable&autowater'};             % right hits, no stim, aw on   (10)
params.condition(end+1) = {'L&hit&~stim.enable&autowater'};             % left hits, no stim, aw on    (11)
% params.condition(end+1) = {'R&miss&~stim.enable&autowater'};            % error right, no stim, aw on  (12)
% params.condition(end+1) = {'L&miss&~stim.enable&autowater'};            % error left, no stim, aw on   (13)

params.tmin = -2.4;
params.tmax = 2.5;
params.dt = 1/50;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect';

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
% % meta = loadEKH1_ALMVideo(meta,datapth); % not usable b/c no usable left miss trials
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB13_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading kinematics ' num2str(sessix) '/' num2str(numel(meta))])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    % kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end



%% Null and Potent Space

clearvars -except obj meta params me sav datapth kin rt

nNPDims = 6;% nNPDims as last arg if you want to change nDims per subspace [4,6,10,13]

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
for sessix = 1:numel(meta)
    disp([meta(sessix).anm ' ' meta(sessix).date])
    % -- input data
    % trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    trialdat_zscored = permute(obj(sessix).trialdat, [1 3 2]);

    % -- null and potent spaces
    cond2use = [1:4 10:11]; % right hit, left hit, right miss, left miss, 2afc
    cond2proj = 1:9; % ~early versions
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0;
    rez(sessix) = singleTrial_elsayed_np_permuteCov(trialdat_zscored, obj(sessix), me(sessix), ...
        params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly); % nNPDims as last arg if you want to change nDims per subspace
   
    % -- coding dimensions
    cond2use = [5 6]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:9]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [5 6]; % for calculating selectivity explained in full neural pop
    rampcond = 9; % corresponding to cond2proj in null/potent analysis
    %     cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,rez(sessix).N_null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
    %     cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,rez(sessix).N_potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);

    cd_null(sessix) = getCodingDimensions(rez(sessix).recon_psth.null,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).recon_psth.potent,trialdat_zscored,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj, rampcond);

end
disp('DONE')

%%

% plotVarianceExplained_main % plot variance explained plots
alignmentPlots_main % plot Alignment plots for NP, Choice, Ramp

%%






























