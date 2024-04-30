clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here

params.region = 'any'; % 'alm','tjm1','mc', 'any'
params.probeType = 'any'; % 'h2','np2','np1', 'any'

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
% meta = loadJPV11(meta,datapth);
meta = loadJPV12(meta,datapth);
% meta = loadJPV13(meta,datapth);

% subset meta
meta = meta(1);
% meta = subsetMetaByParams(meta,params);

%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

me = loadMotionEnergy(obj, meta, params, datapth);

kin = getKinematics(obj, me, params);

% TAGGED UNIT META

tag.nTag = numel(obj.tag);
tag.cluid.clu = [obj.tag(:).cluid]; % where tagged units are in obj.clu
tag.cluid.obj = find(ismember(params.cluid,tag.cluid.clu))'; % where in params.cluid, trialdat, psth

%% SELECTIVITY

% doing this with PSTHs for now, need to do it with single trials (balanced)

%<XY>: X=trial type, Y=lick direction

% conds
cond.RR = 2;
cond.LL = 3;
cond.RL = 4;
cond.LR = 5;

% trials per condition
fnames = fieldnames(cond);
for i = 1:numel(fnames)
    f = fnames{i};
    trials.(f) = params.trialid{cond.(f)};
end

% stimulus selectivity = <RR>+<RL> - <LL>+<LR>
sel.stim = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.RL) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.LR) );

% choice selectivity = <RR>+<LR> - <LL>+<RL>
sel.choice = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LR) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.RL) );

% action selectivity = <RR>+<LR> - <LL>+<RL>
sel.action = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LR) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.RL) );

% outcome selectivity = <LL>+<RR> - <LR>+<RL>
sel.outcome = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LL) ) - ...
    ( obj.psth(:,:,cond.LR)+obj.psth(:,:,cond.RL) );

close all
f = figure;
ax = prettifyAxis(gca);
hold on;
fnames = fieldnames(sel);
cols = linspecer(numel(fnames));
for i = 1:numel(fnames)
    f = fnames{i};
    mu = mean(sel.(f),2);
    ci = getCI(sel.(f));
    shadedErrorBar(obj.time,mu,ci,{'Color',cols(i,:),'LineWidth',2,'HandleVisibility','off'},0.2,ax)
end
xlim([obj.time(1),obj.time(end)])
plotEventTimes(ax,params.eventTimes)


%% find selective neurons
clear times mu sel

%<XY>: X=trial type, Y=lick direction

nboots = 1000;

% conds
cond.RR = 2;
cond.LL = 3;
cond.RL = 4;
cond.LR = 5;

% times (relative to go cue)
times.stim    = [-1.85 -1.2];
times.choice  = [-0.6 0];
times.action  = [0 0.6];
times.outcome = [2 2.995];
fnames = fieldnames(times);
for i = 1:numel(fnames)
    f = fnames{i};
    times.ix.(f) = findTimeIX(obj.time,times.(f), true);
end

for iboot = 1:nboots
    if mod(iboot,50)==0
        disp(['Iteration: ' num2str(iboot) '/' num2str(nboots)])
    end
    % balance # of correct trials and balance # of error trials
    correct_trials = balanceAndSplitTrials(params.trialid,[cond.RR,cond.LL],1.0,0);
    trials.RR = correct_trials.train{1};
    trials.LL = correct_trials.train{2};

    error_trials = balanceAndSplitTrials(params.trialid,[cond.RL,cond.LR],0.9,0.1);
    trials.RL = error_trials.train{1};
    trials.LR = error_trials.train{2};

    % stim
    fnames = fieldnames(trials);
    for icond = 1:numel(fnames)
        f = fnames{icond};
        thesetrials = trials.(f);

        N = numel(thesetrials);

        T = numel(times.ix.stim);

        fr = obj.trialdat(times.ix.stim,:,thesetrials);
        fr_sum_trials = squeeze(sum(fr,3));
        fr_sum_time_trials = sum(fr_sum_trials,1)';

        mu.(f) = 1/(N*T) * fr_sum_time_trials; % (neurons,1)
    end
    
    % sel.stim(:,iboot,1) = mu.RR+mu.RL; % (neurons,boots,group)
    % sel.stim(:,iboot,2) = mu.LL+mu.LR; % (neurons,boots,group)
    
    sel.stim(:,iboot,1) = mu.RR; % (neurons,boots,group)
    sel.stim(:,iboot,2) = mu.LL; % (neurons,boots,group)

    cd(:,iboot) = (mu.RR+mu.RL) - (mu.LL+mu.LR);
    % cd(:,iboot) = mu.RR - mu.LL;

end
for i = 1:numel(params.cluid)
    % [p(i),h(i)] = ranksum(squeeze(sel.stim(i,:,1)) , squeeze(sel.stim(i,:,2)) , 'alpha',0.01 );
    [h(i),p(i)] = ttest(squeeze(sel.stim(i,:,1)) , squeeze(sel.stim(i,:,2)) , 'alpha',0.001);
end
h(isnan(h)) = 0;
find(~h)



%%
a = squeeze(sel.stim(9,:,:));
figure; plot(a)

%%

cd_ = nanmean(cd,2); % mean over boots
% cd_ = cd(:,1);
proj = tensorprod(obj.trialdat,cd_,2,1);


r = proj(:,params.trialid{2});
l = proj(:,params.trialid{3});

mur = mean(r,2);
mul = mean(l,2);
er = getCI(r);
el = getCI(l);

f = figure;
ax = prettifyAxis(gca);
hold on;
shadedErrorBar(obj.time,mur,er,{'Color','b','LineWidth',2},0.2,ax);
shadedErrorBar(obj.time,mul,el,{'Color','r','LineWidth',2},0.2,ax);




