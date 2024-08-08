clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'scripts')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'firstLick';

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
% meta = loadJPV13(meta,datapth);

meta = meta(1);

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
% ------------------------------------------
% -- Kinematics --
% kin (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end

% TAGGED UNIT META

tag.nTag = numel(obj.tag);
tag.cluid = [obj.tag(:).cluid]; % where tagged units are in obj.clu
tag.cluid = find(ismember(params.cluid,tag.cluid))'; % where in params.cluid, trialdat, psth


% % % % % %% CALCULATE CD AT EACH TIME POINT
% % % % %
% % % % % par.train = 0.5;
% % % % % par.test = 1 - par.train;
% % % % % par.cond2use = [2 3]; % rhit lhit
% % % % %
% % % % % [cd.W,cd.data,cd.par] = calcCDAll(obj,params,par);
% % % % %
% % % % % %% PLOT CD
% % % % %
% % % % % cols = getColors;
% % % % % c(1,:) = cols.rhit;
% % % % % c(2,:) = cols.lhit;
% % % % % alph = 0.2;
% % % % %
% % % % % f = figure;
% % % % % f.Renderer = 'painters';
% % % % % ax = prettifyPlot(gca);
% % % % % hold on
% % % % % for i = 1:numel(cd.par.cond2use)
% % % % %     mu = nanmean(cd.data.proj(:,:,i),2);
% % % % %     ci = getCI(cd.data.proj(:,:,i));
% % % % %     shadedErrorBar(obj.time,mu,ci,{'LineWidth',2,'Color',c(i,:)},alph,ax);
% % % % % end
% % % % % plotEventTimes(ax,params.eventTimes);
% % % % % xlim([params.tmin,params.tmax])


%% CALCULATE CODING DIRECTIONS


inp.cond2use = [2 3];

inp.cd_labels = {'stimulus','choice','action'};
inp.cd_epochs = {'goCue','goCue','goCue'};
inp.cd_times = {[-1.85, -1.2],
    [-0.4,0],
    [0 0.4]}; % in seconds, relative to respective epochs

rez = getCodingDimensions_v2(inp,obj,params);

%% plot CDs

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;
alph = 0.2;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
for icd = 1:numel(inp.cd_labels)
    ax = prettifyPlot(nexttile);
    hold on
    for icond = 1:numel(inp.cond2use)
        this = rez.proj{icond}(:,:,icd);
        mu = mean(this,2,'omitnan');
        ci = getCI(this);
        shadedErrorBar(obj.time,mu,ci,{'LineWidth',2,'Color',c(icond,:)},alph,ax)
    end
    plotEventTimes(ax,params.eventTimes);
    xlim([params.tmin,params.tmax])
    title(inp.cd_labels{icd})
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Projection (a.u)')

%% PREDICT NEURAL ACTIVITY FROM KINEMATICS

decodeNeuralFromKin

%% PLOT
clear mu ci
plotDataAndPred

%% create new data objs
% newobj contains trialdat and psth for zscored data
% newobj_move contains trialdat and psth for movement subtracted data

clear sel newobj newobj_move

for i = 1:size(y,3)
    tempy = y(:,:,i);
    tempyhat = yhat(:,:,i);

    % put tempy in same scale as tempyhat
    yhatstd = std(tempyhat,[],1);
    % tempyhat = tempyhat ./ yhatstd;

    % % make sure the baselines start at 0
    % basey = mean(tempy(bix,:),1);
    % tempy = tempy - basey;
    % baseyhat = mean(tempyhat(bix,:),1);
    % tempyhat = tempyhat - baseyhat;

    y_(:,:,i) = tempy;
    yhat_(:,:,i) = tempyhat;
end

newobj = obj;
newobj_move = obj;
newobj.trialdat = permute(y_,[1 3 2]);
newobj_move.trialdat = newobj.trialdat - permute(yhat_, [1 3 2]); % remove predictions
% newobj_move.trialdat = permute(yhat_, [1 3 2]); % remove predictions
for i = 1:numel(params.trialid)
    temp = newobj.trialdat(:,:,params.trialid{i});
    newpsth(:,:,i) = squeeze(mean(temp,3));
    temp = newobj_move.trialdat(:,:,params.trialid{i});
    newpsth_move(:,:,i) = squeeze(mean(temp,3));
end
newobj.psth = newpsth;
newobj_move.psth = newpsth_move;

newobj = rmfield(newobj,{'clu','traj','me','tag','metrics'});
newobj_move = rmfield(newobj_move,{'clu','traj','me','tag','metrics'});
clear newpsth newpsth_move

%% calc selectivity and label cells as sample selective and/or delay selective
clear sel

cond2use = [2 3]; % right hit - left hit
inp.nTrials = 80;
inp.pval = 0.001;
epoch = 'goCue';

whenSelective = [-1.85 -1.2]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.sample.ts,sel.sample.cluid,shadeix.sample,sel.sample.dprime] = ...
    calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [-0.4 -0.0]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.delay.ts,sel.delay.cluid,shadeix.delay,sel.delay.dprime] = ...
    calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [0.0 0.4]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.response.ts,sel.response.cluid,shadeix.response,sel.response.dprime] = ...
    calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

% movement removed data
whenSelective = [-1.85 -1.2]; % seconds relative to epoch, used to find cond a neuron is selective for
[selmove.sample.ts,selmove.sample.cluid,shadeix.sample,selmove.sample.dprime] = ...
    calcPrefSelectivity(newobj_move, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [-0.4 -0.0]; % seconds relative to epoch, used to find cond a neuron is selective for
[selmove.delay.ts,selmove.delay.cluid,shadeix.delay,selmove.delay.dprime] = ...
    calcPrefSelectivity(newobj_move, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [0.0 0.4]; % seconds relative to epoch, used to find cond a neuron is selective for
[selmove.response.ts,selmove.response.cluid,shadeix.response,selmove.response.dprime] = ...
    calcPrefSelectivity(newobj_move, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

%% plot selectivity
close all
alph = 0.15;
cols = linspecer(3);

clear mu
clear ci

% plot selectivity across all units
f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
mu.sample = mean(sel.sample.ts,2);
ci.sample = getCI(sel.sample.ts);
shadedErrorBar(obj.time,mu.sample,ci.sample,{'Color',cols(1,:),'LineWidth',2},alph,ax)
mu.delay = mean(sel.delay.ts,2);
ci.delay = getCI(sel.delay.ts);
shadedErrorBar(obj.time,mu.delay,ci.delay,{'Color',cols(2,:),'LineWidth',2},alph,ax)
mu.response = mean(sel.response.ts,2);
ci.response = getCI(sel.response.ts);
shadedErrorBar(obj.time,mu.response,ci.response,{'Color',cols(3,:),'LineWidth',2},alph,ax)
plotEventTimes(ax,params.eventTimes)
title('DATA')

clear mu
clear ci

% plot selmoveectivity across all units (move removed)
f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
mu.sample = mean(selmove.sample.ts,2);
mu.sample = mu.sample - mean(mu.sample(bix));
ci.sample = getCI(selmove.sample.ts);
shadedErrorBar(obj.time,mu.sample,ci.sample,{'Color',cols(1,:),'LineWidth',2},alph,ax)
mu.delay = mean(selmove.delay.ts,2);
mu.delay = mu.delay - mean(mu.sample(bix));
ci.delay = getCI(selmove.delay.ts);
shadedErrorBar(obj.time,mu.delay,ci.delay,{'Color',cols(2,:),'LineWidth',2},alph,ax)
mu.response = mean(selmove.response.ts,2);
mu.response = mu.response - mean(mu.sample(bix));
ci.response = getCI(selmove.response.ts);
shadedErrorBar(obj.time,mu.response,ci.response,{'Color',cols(3,:),'LineWidth',2},alph,ax)
plotEventTimes(ax,params.eventTimes)
title('MOVE')

%% find how much selectivity per epoch, per unit between move removed and data
close all


epochs = {'sample','delay','response'};
f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
cols = linspecer(3);
c = [252, 248, 3]./255;

allix = 1:numel(params.cluid);
tagix = tag.cluid;
nontagix = allix(~ismember(allix,tagix));

for iepoch = 1:numel(epochs)
    ax = prettifyPlot(nexttile,15);
    hold on;

    shade = shadeix.(epochs{iepoch});

    this = mean(sel.(epochs{iepoch}).ts(shade,nontagix),1);
    thismove = mean(selmove.(epochs{iepoch}).ts(shade,nontagix),1);
    scatter(this,thismove,20,'filled','markeredgecolor','w','linewidth',0.1,'markerfacecolor',cols(iepoch,:))

    this = mean(sel.(epochs{iepoch}).ts(shade,tagix),1);
    thismove = mean(selmove.(epochs{iepoch}).ts(shade,tagix),1);
    scatter(this,thismove,40,'filled','markeredgecolor','k','markerfacecolor',c)

    title(epochs{iepoch})
    mn = min(ax.XLim(1),ax.YLim(1));
    mx = max(ax.XLim(2),ax.YLim(2));
    ax.XLim = [mn,mx];
    ax.YLim = [mn,mx];
    xlims = ax.XLim;
    ylims = ax.YLim;

    ll = line([xlims(1) xlims(2)],[ylims(1) ylims(2)]);
    ll.Color = 'k';
    ll.LineStyle = '--';
    ll.LineWidth = 2;

end

xlabel(t,'Selectivity (zscored)')
ylabel(t,'Selectivity after subtraction (zscored)')


%% CALC CDS FROM NEWOBJ AND NEWOBJ_MOVE

inp.cond2use = [2 3];

inp.cd_labels = {'stimulus','choice','action'};
inp.cd_epochs = {'goCue','goCue','goCue'};
inp.cd_times = {[-1.85, -1.2],
    [-0.4,0],
    [0 0.4]}; % in seconds, relative to respective epochs

rez = getCodingDimensions_v2(inp,newobj,params);

rez_move = getCodingDimensions_v2(inp,newobj_move,params);

%% plot CDs

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;
alph = 0.2;

toplot = rez;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
for icd = 1:numel(inp.cd_labels)
    ax = prettifyPlot(nexttile);
    hold on
    for icond = 1:numel(cond2plot)
        this = toplot.proj{icond}(:,:,icd);
        mu = mean(this,2,'omitnan');
        ci = getCI(this);
        shadedErrorBar(obj.time,mu,ci,{'LineWidth',2,'Color',c(icond,:)},alph,ax)
    end
    plotEventTimes(ax,params.eventTimes);
    xlim([params.tmin,params.tmax])
    title(inp.cd_labels{icd})
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Projection (a.u)')
title(t,'data')

toplot = rez_move;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
for icd = 1:numel(inp.cd_labels)
    ax = prettifyPlot(nexttile);
    hold on
    for icond = 1:numel(cond2plot)
        this = toplot.proj{icond}(:,:,icd);
        mu = mean(this,2,'omitnan');
        ci = getCI(this);
        shadedErrorBar(obj.time,mu,ci,{'LineWidth',2,'Color',c(icond,:)},alph,ax)
    end
    plotEventTimes(ax,params.eventTimes);
    xlim([params.tmin,params.tmax])
    title(inp.cd_labels{icd})
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Projection (a.u)')
title(t,'move')




















