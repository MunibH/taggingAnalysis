clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'facemap')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% % TODO
% - make sure there's consistency between this script and saveDataForFacemap.m
%     - currently, just need to make sure params and meta are set the same

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

params.behav_only = 0;

params.smooth = 3;

params.dt = 1/100;


%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)


%% LOAD DATA

thismeta = meta(1);

[obj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(obj,sesspar,thismeta);

me = loadMotionEnergy(obj,thismeta,params,datapth);
kin = getKinematics(obj,me,params);

%% Visualize ramping


%% Plot all PSTHs in figures
close all

cond2plot = [2 3];
trials = sesspar.trialid(cond2plot);
nTrialsCond = cell2mat(cellfun(@numel,trials(1:end-1),'uni',0));

tix = findTimeIX(obj.time,[sesspar.eventTimes.goCue-0.5 sesspar.eventTimes.goCue],1);

cols = getColors;

xl = [-2.1,params.tmax];

psthsm = 21;
trialsm = 21;
mesm = 11;

feat = 'motion_energy';
ifeat = ismember(kin.featLeg,feat);
medat = kin.dat(:,:,ifeat);

f = figure;
f.Position = [680         484        1059         394];
f.Renderer = 'painters';

ax1 = prettifyAxis(subplot(1,3,1));
pos = get(ax1, 'Position'); % Get current position
pos = [pos(1) pos(2)*1.5 pos(3) pos(4)/2];
set(ax1, 'Position', pos); % Set the new position
hold on;
ax2 = subplot(1,3,2);
hold on;
ax3 = subplot(1,3,3);
hold on;

for iunit = 1:numel(sesspar.cluid{1})
    cla(ax1)
    cla(ax2)
    cla(ax3)
    
    
    
    rdata = squeeze(obj.trialdat{1}(:,iunit,trials{1}));
    rpsth = mean(rdata,2);
    ldata = squeeze(obj.trialdat{1}(:,iunit,trials{2}));
    lpsth = mean(ldata,2);
    
    [~,rix] = sort(mean(rdata(tix,:),1));
    [~,lix] = sort(mean(ldata(tix,:),1));

    trialdata = cat(2,rdata(:,rix),ldata(:,lix));

    rme = medat(:,trials{1});
    lme = medat(:,trials{2});
    mecat = cat(2,rme(:,rix),lme(:,rix));
    
    plot(ax1,obj.time,mySmooth(rpsth,psthsm,'reflect'),'Color',cols.rhit,'LineWidth',2)
    plot(ax1,obj.time,mySmooth(lpsth,psthsm,'reflect'),'Color',cols.lhit,'LineWidth',2)
    plotEventTimes(ax1,sesspar.eventTimes)
    xlim(ax1,[-2.1 obj.time(end)])

    imagesc(ax2,obj.time,1:size(trialdata,2),mySmooth(trialdata,trialsm,'reflect')')
    plotEventTimes(ax2,sesspar.eventTimes,'w')
    for i = 1:numel(nTrialsCond)
        plot(ax2,ax2.XLim,[nTrialsCond(i) nTrialsCond(i)],'w--')
    end
    xlim(ax2,[-2.1 obj.time(end)])
    ylim(ax2,[-0.1 size(trialdata,2)+0.1])
    colormap(ax2,viridis)
    
    
    imagesc(ax3,obj.time,1:size(trialdata,2),mySmooth(mecat,mesm,'reflect')')
    plotEventTimes(ax3,sesspar.eventTimes,'w')
    for i = 1:numel(nTrialsCond)
        plot(ax3,ax3.XLim,[nTrialsCond(i) nTrialsCond(i)],'w--')
    end
    xlim(ax3,[-2.1 obj.time(end)])
    ylim(ax3,[-0.1 size(trialdata,2)+0.1])
    % colormap(ax3,flipud(rocket))
    % colorbar
    % clim([200 180])

break
end
% xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',18)
% ylabel(t,'Spks/sec','FontSize',18)