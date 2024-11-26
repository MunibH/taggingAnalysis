clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'medullaDynamics')));

clc

%% PARAMETERS

% load data params
params = defaultParams();
params.dt = 1/100;
params.behav_only = 1;

% model params
p.hidden_dim = 128;
p.num_layers = 3;
p.Animal = 0;
p.model_params_filepath = fullfile(utilspth,'medullaDynamics','trained_model_parameters_20241112.pth');
p.nPCs = 8; % eval('size(kin.dat,3)');
p.overwrite = false; % if 1 and results already exist, overwrite. if 0, skip if results already exist

%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
savepth = fullfile(datapth,'ImputedMedullaDynamics');
if ~exist(savepth)
    mkdir(savepth)
end
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

%% IMPUTE MEDULLA DYNAMICS USING LSTM MODEL

clear rez
for isess = 1:numel(meta)
    clearvars -except meta isess p params datapth savepth utilspth rez

    % set save path 
    p.savepth = fullfile(savepth,[meta(isess).anm '_' meta(isess).date]);
    
    % check if results exist
    cont = 1;
    if p.overwrite
        cont = 1;
    elseif ~p.overwrite && exist([p.savepth '_PCs.csv'], 'file')
        cont = 0;
    elseif ~p.overwrite && ~exist([p.savepth '_PCs.csv'], 'file')
        cont = 1;
    end

    if cont
        % load data
        [sessobj,sesspar] = loadSessionData(meta(isess),params);
        me = loadMotionEnergy(sessobj, meta(isess), sesspar, datapth);
        kin = getKinematics(sessobj, me, sesspar);
        
        % get trial length
        p.trial_len = size(kin.dat,1);
        p.num_trials = size(kin.dat,2);
    
        % reshape kinematic data and save to csv file 
        temp = reshape(kin.dat,p.trial_len*p.num_trials,size(kin.dat,3));
        p.kinpth = [p.savepth '_kinematics.csv'];
        csvwrite(p.kinpth, temp);
    
        medPCs = PyCallLSTM(p); % (time,trials,pcs,session)
    else
        [sessobj,sesspar] = loadSessionData(meta(isess),params);
        me = loadMotionEnergy(sessobj, meta(isess), sesspar, datapth);
        kin = getKinematics(sessobj, me, sesspar);
        pc = csvread([p.savepth '_PCs.csv']);  
        medPCs = reshape(pc, numel(sessobj.time), [], size(pc,2)); % (time,trials,pcs,session)
    end

    
    % ~~~~~~~~~~~PLOT FIRST PC AND A KINEMATIC FEATURE~~~~~~~~~~~
    featstr = 'jaw_ydisp_view1';
    % featstr = 'motion_energy';
    % PlotJawAndPC(medPCs(:,:,1), ...
    %     kin.dat(:,:,ismember(kin.featLeg,featstr)), ...
    %     meta(isess), sessobj.time, sesspar.eventTimes, featstr)

    
    % ~~~~~~~~~~~PLOT SAME THING BY CONDITION~~~~~~~~~~~
    trials = findTrials(sessobj, {'R&hit&~autowater&~early', 'L&hit&~autowater&~early'});
    featstr = 'jaw_ydisp_view1';
    % featstr = 'motion_energy';
    % PlotJawAndPC(medPCs(:,:,1), ...
    %     kin.dat(:,:,ismember(kin.featLeg,featstr)), ...
    %     meta(isess), sessobj.time, sesspar.eventTimes, featstr, trials)

    % ~~~~~~~~~~~SELECTIVITY~~~~~~~~~~~
    trials = findTrials(sessobj, {'R&hit&~autowater&~early', 'L&hit&~autowater&~early'});
    featstr = 'jaw_ydisp_view1';
    % featstr = 'motion_energy';
    % PlotJawAndPC_Selectivity(medPCs(:,:,1), ...
    %     kin.dat(:,:,ismember(kin.featLeg,featstr)), ...
    %     meta(isess), sessobj.time, sesspar.eventTimes, featstr, trials)

    % ~~~~~~~~~~~FEATURE CORRELATION~~~~~~~~~~~
    feats = {'jaw_ydisp_view1','motion_energy'};
    % featstr = 'motion_energy';
    rez.cc(:,:,isess) = CorrelateFeatAndPC(medPCs, kin.dat(:,:,ismember(kin.featLeg,feats))); % (pc,feat,session)

    % ~~~~~~~~~~~FEATURE CROSS-CORRELATION~~~~~~~~~~~
    feats = {'jaw_ydisp_view1','motion_energy'};
    maxLag = 30; % bins
    [rez.xc(:,:,:,isess), lags, peakLags] = CrossCorrelateFeatAndPC(medPCs, kin.dat(:,:,ismember(kin.featLeg,feats)), maxLag); % (lag,pc,feat,session), lags, (pc,feat,session)
    rez.peakLags(:,:,isess) = lags(peakLags).*params.dt.*1000; % in ms, (pc,feat,session)

end

%% REZ

close all

% PLOT CORR

feats = {'Jaw disp','Motion energy'};

temp = permute(rez.cc,[2 1 3]); % (feat,pc,session)

xs = 1:size(temp,2);
c = internet(8);
for ifeat = 1:size(temp,1)
    f = figure;
    f.Position = [680   612   342   266];
    ax = prettifyAxis(gca);
    hold on;
    for ipc = 1:size(temp,2)
        this = squeeze(temp(ifeat,ipc,:));
        x = simple_violin_scatter(xs(ipc)*ones(size(this)), this, numel(this)./ipc, 0.5);
        scatter(x, this, 20,'filled', 'markerfacecolor',c(ipc,:), 'markeredgecolor','none')
    end
    xlabel('PC')
    ylabel('Correlation')
    title(feats{ifeat})
end

% PLOT XCORR
temp = permute(rez.xc,[1 3 2 4]); % (lag,feat,pc,session)
xs = 1:size(temp,2);
c = internet(8);
for ifeat = 1:size(temp,2)
    f = figure;
    f.Position = [680   612   342   266];
    ax = prettifyAxis(gca);
    hold on;
    for ipc = 1:size(temp,3)
        this = squeeze(temp(:,ifeat,ipc,:));
        mu = nanmean(this,2);
        [mu,sig]= mean_CI(this,0.95);
        shadedErrorBar(lags.*params.dt * 1000,mu,sig,{'Color',c(ipc,:),'LineStyle','-','LineWidth',1.5},0.3,ax)
        % plot(lags.*params.dt * 1000,this,'Color',c(ipc,:),'LineStyle','-','LineWidth',1.5)
    end
    xlabel('Lag (ms)')
    ylabel('Correlation')
    title(feats{ifeat})
    xlim([-300 300])
end

% PLOT PEAK XCORR LAG
temp = permute(rez.peakLags,[2 1 3]); % (feat,pc,session)
xs = 1:size(temp,2);
c = internet(8);
for ifeat = 1:size(temp,1)
    f = figure;
    f.Position = [680   612   342   266];
    ax = prettifyAxis(gca);
    hold on;
    for ipc = 1:size(temp,2)
        this = squeeze(temp(ifeat,ipc,:));
        x = simple_violin_scatter(xs(ipc)*ones(size(this)), this, numel(this)./ipc, 0.5);
        scatter(x, this, 20,'filled', 'markerfacecolor',c(ipc,:), 'markeredgecolor','none')
    end
    xlabel('PC')
    ylabel('Peak XCorr Lag (ms)')
    title(feats{ifeat})
end


%% Helper functions

%% plot first PC for all sessions
function PlotJawAndPC(firstpc, jaw, meta, time, evtimes, featstr, varargin)

if nargin>6
    trials = varargin{1};
    alltrials = cell2mat(trials');
    trialix = 1:numel(alltrials);
else
    alltrials = 1:size(jaw,2);
    trialix = alltrials;
end

f = figure;
f.Position = [680   480   687   398];
f.Renderer = 'painters';

ax = prettifyAxis(subplot(1,2,1));
hold on;
this = jaw(:,alltrials);
% this = removeOutliers(this,2.5);
this = mySmooth(this,11,'reflect');
imagesc(time,trialix,this');
clim([135 180])
ax.YDir = 'reverse';
ax.YLim = [0 trialix(end)];
ax.XLim = [time(10) time(end)];
plotEventTimes(ax,evtimes,'w')
colormap(ax,viridis)
title(ax,strrep(featstr,'_',' '));
% xlabel('Time from go cue (s)')
ylabel('Trials')
if nargin>6
    ctrials = cumsum(cell2mat(cellfun(@numel,trials,'uni',0)));
    for i = 1:(numel(trials)-1)
        plot(ax.XLim,[ctrials(i) ctrials(i)], 'w-')
    end
end


ax = prettifyAxis(subplot(1,2,2));
hold on;
this = firstpc(:,alltrials);
imagesc(time,trialix,this');
ax.YDir = 'reverse';
ax.YLim = [0 trialix(end)];
ax.XLim = [time(10) time(end)];
plotEventTimes(ax,evtimes,'w')
colormap(ax,magma)
title(ax,'LSTM PC1')
if nargin>6
    ctrials = cumsum(cell2mat(cellfun(@numel,trials,'uni',0)));
    for i = 1:(numel(trials)-1)
        plot(ax.XLim,[ctrials(i) ctrials(i)], 'w-')
    end
end


sgtitle([meta.anm ' ' meta.date])

drawnow;

end


%% selectivity

function PlotJawAndPC_Selectivity(firstpc, jaw, meta, time, evtimes, featstr, trials)

alltrials = cell2mat(trials');
trialix = 1:numel(alltrials);

cols = getColors();

f = figure;
f.Position = [680   603   779   275];
f.Renderer = 'painters';

ax = prettifyAxis(subplot(1,2,1));
hold on;
this = cell2mat(cellfun(@(x)nanmean(jaw(:,x),2),trials,'uni',0));
plot(time,this(:,1),'Color',cols.rhit,'LineWidth',1.5);
plot(time,this(:,2),'Color',cols.lhit,'LineWidth',1.5);
ax.XLim = [time(10) time(end)];
plotEventTimes(ax,evtimes,'k')
% xlabel('Time from go cue (s)')
ylabel(strrep(featstr,'_',' '))



ax = prettifyAxis(subplot(1,2,2));
hold on;
this = cell2mat(cellfun(@(x) nanmean(firstpc(:,x),2),trials,'uni',0));
plot(time,this(:,1),'Color',cols.rhit,'LineWidth',1.5);
plot(time,this(:,2),'Color',cols.lhit,'LineWidth',1.5);
ax.XLim = [time(10) time(end)];
plotEventTimes(ax,evtimes,'k')
% xlabel('Time from go cue (s)')
ylabel('LSTM PC1');

sgtitle([meta.anm ' ' meta.date])

drawnow;

end

%% correlation

function cc = CorrelateFeatAndPC(pc, kin)

x = reshape(pc,[],size(pc,3));
k = reshape(kin,[],size(kin,3));

for ipc = 1:size(x,2)
    for ifeat = 1:size(k,2)
        temp = corrcoef(x(:,ipc),k(:,ifeat));
        cc(ipc,ifeat) = temp(1,2);
    end
end

end

%% xcorr
function [cc, lags, peakLag] = CrossCorrelateFeatAndPC(pc, kin, maxLag)

% % Set the range of lags
% maxLag = 10;
lags = -maxLag:maxLag;

% Reshape the input matrices to 2D for easier manipulation
x = reshape(pc, [], size(pc, 3));  % [time, pcs]
k = reshape(kin, [], size(kin, 3)); % [time, feats]

% Initialize the cross-correlation matrix
cc = zeros(length(lags), size(x, 2), size(k, 2)); % [lags, pcs, feats]

% Loop over each column in pc and kin to compute cross-correlations
for ipc = 1:size(x, 2)
    for ifeat = 1:size(k, 2)
        temp = xcorr(x(:, ipc) - nanmean(x(:,ipc)), k(:, ifeat) - nanmean(k(:,ifeat)), maxLag, 'coeff'); % Cross-correlation for -10 to +10 lags
        cc(:, ipc, ifeat) = temp; % Store the cross-correlation values for each lag
        [~,peakLag(ipc,ifeat)] = max(cc(:,ipc,ifeat));
    end
end

end



%% Call LSTM Model
function [Predicted_PCs,ve] = PyCallLSTM(p)
% p is model parameter struct

% Construct the command to call the batch file with arguments
command = sprintf('RunModel.bat "%s" %d %d %d %d "%s" "%s" "%d" "%d"', ...
    p.savepth, p.hidden_dim, p.num_layers, p.Animal, p.trial_len, p.model_params_filepath, p.kinpth, p.nPCs, p.num_trials);

% Execute the command
[status, result] = system(command);

% Check for errors
if status == 0
    disp('Saved Imputed Medulla PCs');
    disp(result);

    % Load the predicted PCs
    Predicted_PCs = csvread([p.savepth '_PCs.csv']);  
    ve = csvread([p.savepth '_vExp.csv']);  

    % Reshape if needed (based on trial length and number of PCs)
    Predicted_PCs = reshape(Predicted_PCs, p.trial_len, p.num_trials, size(Predicted_PCs,2)); % [time, trials, PCs]
else
    disp('Error in executing Python script');
    disp(result);
    disp(p)
end

end











