clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'facemap')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

params.behav_only = 0;

params.smooth = 51;

params.dt = 1/75;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
resultspth = fullfile(utilspth,'results','FacemapCodingDirection');
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
clear,clc,close all

resultsdate = '20241010';

fetchpth = fullfile('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\FacemapCodingDirection',resultsdate);
contents = dir(fetchpth);
contents = contents(3:end); % remove . ..


plt.sel = nan(numel(contents,3));
for i = 1:numel(contents) % sessions
    disp(['Session: ' num2str(i) '/' num2str(numel(contents))])
    load(fullfile(fetchpth,contents(i).name))
    
    right_trials = sesspar.trialid{2};
    left_trials = sesspar.trialid{3};
    trialtypes = [-1*ones(numel(right_trials),1) ; ones(numel(left_trials),1)];

    X = cat(2, trialproj.data(:,right_trials), trialproj.data(:,left_trials));
    [accuracy.data(:,i), ~] = FacemapDecodeTrialType(trialtypes, X, 10);

    X = cat(2, trialproj.resid(:,right_trials), trialproj.resid(:,left_trials));
    [accuracy.resid(:,i), ~] = FacemapDecodeTrialType(trialtypes, X, 10);

 end

%%
close all
time = params.tmin:params.dt:params.tmax;
time = time(1:end-1);

f = figure;
f.Position = [680   589   443   289];
f.Renderer = 'painters';
ax = prettifyAxis(gca,'def',15);
hold on;

cc = {[88, 52, 235]./255,[0,0,0]};

fns = fieldnames(accuracy);
for i = 1:numel(fns)
    f = fns{i};
    this = accuracy.(f);
    mu = nanmean(this,2);
    ci = mean_CI(this)./size(this,2);
    mu = mu - nanmean(mu(1:12,:));
    shadedErrorBar(time,mu,ci,{'Color',cc{i},'LineStyle','-','LineWidth',2},0.2,ax)
end
xlabel('Time from go cue')
ylabel('Accuracy')
title('Trial type decoding','FontSize',12)
xlim([-2 3])
ll=line([-1.85 -1.85],ax.YLim);ll.LineStyle = '--';
ll=line([-1.2 -1.2],ax.YLim);ll.LineStyle = '--';
ll=line([0 0],ax.YLim);ll.LineStyle = '--';

lh1 = plot([nan nan],[nan nan],'Color',cc{1},'LineWidth',2);
lh2 = plot([nan nan],[nan nan],'Color',cc{2},'LineWidth',2);
% Manually create the legend with only the last two lines
leg = legend([lh1 lh2], 'Data', 'Facemap residuals');
leg.Box = 'off';
leg.FontSize = 10;

