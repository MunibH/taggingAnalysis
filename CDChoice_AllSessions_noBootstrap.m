clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - time warp
% - fix jaw phase
%   - it now is a cell array, need to interp the time like ME


% for early lick trial formatingg
% - the video frame times for 1 or 2 sessions are a little wonky. see n1
% and n2 below.
%               'JPV11'
%          date: '2023-06-21'
% above session, side cam for all trials starts at a time that is greater
% than 0.....
% SOLUTION: figure out if those frame times are delayed w.r.t to trial
% start or if they just dropped frames from trial start to then
% then, write a function that identifies this and aligns neural data to
% video data properly. don't hard code that 0.5 anymore that could mean...

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single'};

params.behav_only = 0;

params.useEarly = 0;

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

%% LOAD DATA

clear cd_corr
for isess = 1:numel(meta)
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);

    cdix = findTimeIX(sessobj.time,[sesspar.eventTimes.goCue-0.6 sesspar.eventTimes.goCue],1); % choice
    normix = findTimeIX(sessobj.time,[-0.4,0],1);
    baselineix = findTimeIX(sessobj.time,[-2 2]);

    % get balanced number of trials
    condix = find(ismember(sesspar.condLabel,{'rhit','lhit'}));
    trials = balanceAndSplitTrials(sesspar.trialid, condix, 1, 0); % {rhit,lhit}
    trials = cellfun(@(x)  randsample(x,50,true), trials, 'uni', 0);
    % trials = sesspar.trialid(condix);

    % sample single units
    trialdat_un = cat(2,sessobj.trialdat{:});
    unitix = randsample(size(trialdat_un,2),50,true);
    trialdat_un = trialdat_un(:,unitix,:);

    % calculate PSTH
    for icond = 1:numel(trials)
        psth(:,:,icond) = nanmean(trialdat_un(:,:,trials{icond}),3);
    end

    % mean subtract PSTH using baseline mean
    mu_un = nanmean(psth(baselineix,:,:),[1 3]);
    mu_un = fillmissing(mu_un,"constant",nanmedian(mu_un));
    mu_un(mu_un==0) = 1;
    sigma_un = nanstd(psth(baselineix,:,:),[],[1 3]);
    sigma_un = fillmissing(sigma_un,"constant",nanmedian(sigma_un));
    sigma_un(sigma_un==0) = 1;

    psth = (psth - mu_un) ./ sigma_un;

    clear pref_dir temp
    for iunit = 1:size(psth,2)
        temp = squeeze(psth(cdix,iunit,:));
        a = temp(:,1) > temp(:,2);
        a = a(:);
        a = sum(a) / numel(a);
        if a < 0.5
            pref_dir(iunit,:) = [2 1];
        else
            pref_dir(iunit,:) = [1 2];
        end
    end
    temp = psth;
    for iunit = 1:size(psth,2)
        temp(:,iunit,1) = psth(:,iunit,pref_dir(iunit,1));
        temp(:,iunit,2) = psth(:,iunit,pref_dir(iunit,2));
    end
    psth = temp;

    % CALCULATE CD CHOICE
    cd = calcCD(psth,cdix);


    % CALCULATE CD EACH TIME POINT
    for itime = 1:size(psth,1)    
        cd_time = calcCD(psth(itime,:,:),1);
        cd_corr(itime,isess) = corr(cd_time,cd);
    end
    
    % cd_corr(:,isess) = cd_corr(:,isess) ./ nanmean(cd_corr(:,isess));
end

%% PLOT

cols = getColors;

f = figure;
f.Position = [680   564   512   314];
ax = prettifyAxis(gca);
hold on;
[m,h] = mean_CI(cd_corr,0.95,false);
shadedErrorBar(sessobj.time*1000,m+0.6,h,{'Color',cols.un,'LineWidth',2},0.3,ax)
xlim([-80 100])
xlabel('Time from go cue (ms)')
ylabel({'Normalized correlation','with CD choice'})
ylim([0 ax.YLim(2)])
plot([0 0],ax.YLim,'k--')













