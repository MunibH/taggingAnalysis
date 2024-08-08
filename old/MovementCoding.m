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
% params.alignEvent = 'firstLick';

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
% meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
meta = loadJPV13(meta,datapth);

meta = meta(2);

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

% event times
params.eventTimes = getEventTimes(obj.bp.ev,params.events,params.alignEvent);


% TAGGED UNIT META

tag.nTag = numel(obj.tag);
tag.cluid = [obj.tag(:).cluid]; % where tagged units are in obj.clu
tag.cluid = find(ismember(params.cluid,tag.cluid))'; % where in params.cluid, trialdat, psth

% GET EVENT TIMES
evtimes = getEventTimes(obj(1).bp.ev,params(1).events,params(1).alignEvent);
align = mode(obj(1).bp.ev.(params(1).alignEvent));

%% calc selectivity and label cells as sample selective and/or delay selective
clear sel

cond2use = [2 3]; % right hit - left hit
inp.nTrials = 50;
inp.pval = 0.001;
epoch = 'goCue';

whenSelective = [-1.85 -1.2]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.sample.ts,sel.sample.cluid,shadeix.sample,sel.sample.dprime] = ...
    calcPrefSelectivity(obj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [-0.4 -0.0]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.delay.ts,sel.delay.cluid,shadeix.delay,sel.delay.dprime] = ...
    calcPrefSelectivity(obj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [0.0 0.4]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.response.ts,sel.response.cluid,shadeix.delay,sel.response.dprime] = ...
    calcPrefSelectivity(obj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

dprime = calcdprime(obj,params,cond2use,inp);

%% plot selectivity
close all
alph = 0.15;
cols = linspecer(3);

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
title('All units')

clear mu
clear ci

% plot selectivity across selective units
f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
temp = sel.sample.ts(:,sel.sample.cluid);
mu.sample = mean(temp,2);
ci.sample = getCI(temp);
shadedErrorBar(obj.time,mu.sample,ci.sample,{'Color',cols(1,:),'LineWidth',2},alph,ax)
temp = sel.delay.ts(:,sel.delay.cluid);
mu.delay = mean(temp,2);
ci.delay = getCI(temp);
shadedErrorBar(obj.time,mu.delay,ci.delay,{'Color',cols(2,:),'LineWidth',2},alph,ax)
temp = sel.response.ts(:,sel.response.cluid);
mu.response = mean(temp,2);
ci.response = getCI(temp);
shadedErrorBar(obj.time,mu.response,ci.response,{'Color',cols(3,:),'LineWidth',2},alph,ax)
plotEventTimes(ax,params.eventTimes)
title('Selective units')

clear mu
clear ci

% plot selectivity across non-selective units
f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
temp = sel.sample.ts(:,~sel.sample.cluid);
mu.sample = mean(temp,2);
ci.sample = getCI(temp);
shadedErrorBar(obj.time,mu.sample,ci.sample,{'Color',cols(1,:),'LineWidth',2},alph,ax)
temp = sel.delay.ts(:,~sel.delay.cluid);
mu.delay = mean(temp,2);
ci.delay = getCI(temp);
shadedErrorBar(obj.time,mu.delay,ci.delay,{'Color',cols(2,:),'LineWidth',2},alph,ax)
temp = sel.response.ts(:,~sel.response.cluid);
mu.response = mean(temp,2);
ci.response = getCI(temp);
shadedErrorBar(obj.time,mu.response,ci.response,{'Color',cols(3,:),'LineWidth',2},alph,ax)
plotEventTimes(ax,params.eventTimes)
title('Non-selective units')

clear mu
clear ci

%% PLOT DPRIME

close all

baselineedges = [params.eventTimes.bitStart+0.1 params.eventTimes.sample];
baselineix = findTimeIX(obj.time,baselineedges);
bix = baselineix(1):baselineix(2);

cols = getColors;
cmap = flipud(redblue);
alph = 0.2;
dprimethresh = 0.0;

dprime_ = dprime * 3;

mn = min(min(dprime_));
mx = max(max(dprime_));

lab = {'sample','delay','response'};
edges = -1:0.1:1;
f = figure; f.Position = [391         227        1199         678];
f.Renderer = 'painters';
% sort by sample epoch dprime
ax = prettifyPlot(subplot(2,3,1));
hold on;
ix = [-1.85 -1.2]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_,2),dprime_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / 5;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:37;
ylim([notSelective(end)+1 size(dprime_,2)])
title('sample')
ylabel('Units')
ax = prettifyPlot(subplot(2,3,4));
hold on;
temp = dprime_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)
ylabel('d-prime')
xlabel('Time from go cue (s)')

% sort by delay epoch dprime
ax = prettifyPlot(subplot(2,3,2));
hold on;
ix = [-0.4 0]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_,2),dprime_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / 5;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:31;
ylim([notSelective(end)+1 size(dprime_,2)])
title('delay')
ax = prettifyPlot(subplot(2,3,5));
hold on;
temp = dprime_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < -dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)

% sort by response epoch dprime
ax = prettifyPlot(subplot(2,3,3));
hold on;
ix = [0 0.4]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_,2),dprime_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / 5;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:69;
ylim([notSelective(end)+1 size(dprime_,2)])
title('response')
ax = prettifyPlot(subplot(2,3,6));
hold on;
temp = dprime_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < -dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)

%% use kinematics to decode neural activity, subtract off portion of neural activity that is predicted by movements
clear par
% DECODING PARAMETERS -----------------------------------------------------

% input data = kin data (time*trials,kinfeats)
% output data = neural data   (time*trials,units)

par.pre=15; % time bins prior to output used for decoding
par.post=2; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% data sets
par.train = 0.7; % fraction of trials
par.test = 1 - par.train;

% feature to use to decode
par.feats = kin(1).featLeg;
par.feats = {'motion','nose','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});

% trials
par.cond2use = [2, 3];
par.conds = {'rhit','lhit'};

% DECODING -----------------------------------------------------

% trials
for i = 1:numel(par.cond2use)
    par.trials.(par.conds{i}) = params.trialid{par.cond2use(i)};
    nTrialsCond(i) = numel(par.trials.(par.conds{i}));
end
% balance number of trials across conditions
minTrials = min(nTrialsCond);
mask = nTrialsCond > minTrials;
par.trials.all = [];
for i = 1:numel(par.cond2use)
    if mask(i)
        par.trials.(par.conds{i}) = sort(randsample(par.trials.(par.conds{i}), minTrials, false),'ascend');
    end
    par.trials.all = [par.trials.all; par.trials.(par.conds{i})];
end

% partition train and test
nTrials = numel(par.trials.all);
nTrain = floor(nTrials*par.train);
par.trials.train = randsample(par.trials.all,nTrain,false);
par.trials.test = par.trials.all(~ismember(par.trials.all,par.trials.train));

% input data
par.featix = find(ismember(kin.featLeg,par.feats));

X.train = kin.dat(:,par.trials.train,par.featix); % (time,trials,feats)
X.size.train = size(X.train);
X.train = reshape(X.train, size(X.train,1)*size(X.train,2),size(X.train,3));

X.test = kin.dat(:,par.trials.test,par.featix); % (time,trials,feats)
X.size.test = size(X.test);
X.test = reshape(X.test, size(X.test,1)*size(X.test,2),size(X.test,3));

% reshape train and test data to account for prediction bin size
X.train = reshapePredictors(X.train,par);
X.test = reshapePredictors(X.test,par);

% flatten inputs
% if you're using a model with recurrence, don't flatten
X.train = reshape(X.train,size(X.train,1),size(X.train,2)*size(X.train,3));
X.test = reshape(X.test,size(X.test,1),size(X.test,2)*size(X.test,3));

% output data
Y.train = permute(obj.trialdat(:,:,par.trials.train), [1 3 2]); % (time,trials,units);
Y.size.train = size(Y.train);
Y.train = reshape(Y.train, size(Y.train,1)*size(Y.train,2),size(Y.train,3));

Y.test = permute(obj.trialdat(:,:,par.trials.test), [1 3 2]); % (time,trials,units);
Y.size.test = size(Y.test);
Y.test = reshape(Y.test, size(Y.test,1)*size(Y.test,2),size(Y.test,3));

% standardize data
% standardize both train and test sets using train set statistics
% can also standardize using specific time points (presample for example)
X.mu = mean(X.train,1,'omitnan');
X.sigma = std(X.train,[],1,'omitnan');
X.train = (X.train - X.mu) ./ X.sigma;
if ~par.test==0
    X.test = (X.test - X.mu) ./ X.sigma;
end

Y.mu = mean(Y.train,1,'omitnan');
Y.sigma = std(Y.train,[],1,'omitnan');
Y.train = (Y.train - Y.mu) ./ Y.sigma;
if ~par.test==0
    Y.test = (Y.test - Y.mu) ./ Y.sigma;
end

% fill missing values in kinematics
X.train = fillmissing(X.train,'constant',0);
Y.train = fillmissing(Y.train,'nearest');
X.test = fillmissing(X.test,'constant',0);
Y.test = fillmissing(Y.test,'nearest');

% train and test models
for i = 1:Y.size.train(3)
    if mod(i,20) == 0
        disp(['Unit ' num2str(i) '/' num2str(Y.size.train(3))])
    end
    mdl{i} = fitrlinear(X.train,Y.train(:,i));
    Y.pred(:,i) = predict(mdl{i},X.test);
end

% variance explained
y = reshape(Y.test,Y.size.test(1),Y.size.test(2),[]); % original input data (centered)
yhat = reshape(Y.pred,Y.size.test(1),Y.size.test(2),[]); % prediction
for i = 1:size(y,3)
    tempcorr = corrcoef(y(:,i),yhat(:,i));
    Y.R2(i) = tempcorr(1,2).^2;
end

%% predict activity for all trials
clear XX YY

XX = kin.dat(:,:,par.featix); % (time,trials,feats)
XX = reshape(XX, size(XX,1)*size(XX,2),size(XX,3));

% reshape train and test data to account for prediction bin size
XX = reshapePredictors(XX,par);

% flatten inputs
% if you're using a model with recurrence, don't flatten
XX = reshape(XX,size(XX,1),size(XX,2)*size(XX,3));

% output data
YY.data = permute(obj.trialdat, [1 3 2]); % (time,trials,units);
YY.data = reshape(YY.data, size(YY.data,1)*size(YY.data,2),size(YY.data,3));

% standardize data
% standardize both train and test sets using train set statistics
% can also standardize using specific time points (presample for example)
XXmu = mean(XX,1,'omitnan');
XXsigma = std(XX,[],1,'omitnan');
XX = (XX - XXmu) ./ XXsigma;

YY.mu = mean(YY.data,1,'omitnan');
YY.sigma = std(YY.data,[],1,'omitnan');
YY.data = (YY.data - YY.mu) ./ YY.sigma;

% fill missing values in kinematics
XX = fillmissing(XX,'constant',0);

% predictions
for i = 1:numel(mdl)
    if mod(i,20) == 0
        disp(['Unit ' num2str(i) '/' num2str(Y.size.train(3))])
    end
    YY.pred(:,i) = predict(mdl{i},XX);
end

y = reshape(YY.data,numel(obj.time),obj.bp.Ntrials,[]); % original input data (centered)
yhat = reshape(YY.pred,numel(obj.time),obj.bp.Ntrials,[]); % prediction (time,trials,units)

% variance explained
for i = 1:size(y,3)
    tempcorr = corrcoef(y(:,i),yhat(:,i));
    YY.R2(i) = tempcorr(1,2).^2;
end

%% make susu fig 3 plots

close all

cols = getColors();

cond2plot = [2,3];
trix = cell2mat(params.trialid(cond2plot)');
for i = 1:numel(cond2plot)
    ctrix{i} = params.trialid{cond2plot(i)};
end

f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
histogram(YY.R2,20,'Normalization','count')
xlabel('Test-R2')
ylabel('Unit count')


baselineedges = [params.eventTimes.bitStart+0.1 params.eventTimes.sample];
baselineix = findTimeIX(obj.time,baselineedges);
bix = baselineix(1):baselineix(2);


f = figure;
f.Renderer = 'painters';
ax1 = subplot(2,2,1);
hold on;
ax2 = subplot(2,2,2);
hold on;
ax3 = subplot(2,2,3);
hold on;
ax4 = subplot(2,2,4);
hold on;
for i = 1:size(y,3)
    cla(ax1)
    cla(ax2)
    cla(ax3)
    cla(ax4)

    % if ~sel.delay.cluid(i)
    %     continue
    % end

    tempy = y(:,:,i);
    tempyhat = yhat(:,:,i);

    % put tempy in same scale as tempyhat
    yhatstd = std(tempyhat,[],1);
    tempyhat = tempyhat ./ yhatstd;
    tempy = tempy;% .* yhatstd;

    % make sure the baselines start at 0
    basey = mean(tempy(bix,:),1);
    tempy = tempy - basey;
    baseyhat = mean(tempyhat(bix,:),1);
    tempyhat = tempyhat - baseyhat;


    hold(ax1,'on')
    imagesc(ax1,obj.time,1:size(y,2),tempy');
    ax1.YDir = "reverse";
    title(ax1,'data')
    plotEventTimes(ax1,params.eventTimes)
    ylim(ax1,[1 size(y,2)])
    colorbar(ax1);

    hold(ax2,'on')
    imagesc(ax2,obj.time,1:size(y,2),tempyhat');
    ax2.YDir = "reverse";
    title(ax2,'pred')
    plotEventTimes(ax2,params.eventTimes)
    ylim(ax2,[1 size(y,2)])
    colorbar(ax2);

    
    tempr = tempy(:,ctrix{1});
    templ = tempy(:,ctrix{2});
    temprhat = tempyhat(:,ctrix{1});
    templhat = tempyhat(:,ctrix{2});

    mu.y{1} = mean(tempr,2);
    mu.y{2} = mean(templ,2);
    mu.yhat{1} = mean(temprhat,2);
    mu.yhat{2} = mean(templhat,2);

    ci.y{1} = getCI(tempr);
    ci.y{2} = getCI(templ);
    ci.yhat{1} = getCI(tempr);
    ci.yhat{2} = getCI(templ);

    hold(ax3,'on')
    shadedErrorBar(obj.time,mu.y{1},ci.y{1},{'Color',cols.rhit,'LineWidth',2},alph,ax3)
    shadedErrorBar(obj.time,mu.y{2},ci.y{2},{'Color',cols.lhit,'LineWidth',2},alph,ax3)
    hold(ax4,'on')
    shadedErrorBar(obj.time,mu.yhat{1},ci.yhat{1},{'Color',cols.rhit_aw,'LineWidth',2},alph,ax4)
    shadedErrorBar(obj.time,mu.yhat{2},ci.yhat{2},{'Color',cols.lhit_aw,'LineWidth',2},alph,ax4)
    title(ax4,['R2=' num2str(round(YY.R2(i),2))])
    

    % ylims(1) = min(ax3.YLim(1),ax4.YLim(1));
    % ylims(2) = max(ax3.YLim(2),ax4.YLim(2));
    % % ylims(1) = ylims(1)*1.1;
    % % ylims(2) = ylims(2)*1.1;
    % 
    % ax3.YLim = ylims;
    % ax4.YLim = ylims;

    plotEventTimes(ax3,params.eventTimes)
    plotEventTimes(ax4,params.eventTimes)
    xlabel(ax3,'Time from go cue (s)')
    ylabel(ax3,'zscored FR')
    xlabel(ax4,'Time from go cue (s)')
    ylabel(ax4,'zscored FR')
    title(ax3,[meta.anm ' ' meta.date ' Unit ' num2str(i)])

    pause

    clear mu ci
end

% clear mu ci

%% create new data objs
% newobj contains trialdat and psth for zscored data
% newobj_move contains trialdat and psth for movement subtracted data

clear sel newobj newobj_move

for i = 1:size(y,3)
    tempy = y(:,:,i);
    tempyhat = yhat(:,:,i);

    % put tempy in same scale as tempyhat
    yhatstd = std(tempyhat,[],1);
    tempyhat = tempyhat ./ yhatstd;

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
newobj_move.trialdat = y_ - permute(yhat_, [1 3 2]); % remove predictions
for i = 1:numel(params.trialid)
    temp = newobj.trialdat(:,:,params.trialid{i});
    newpsth(:,:,i) = squeeze(mean(temp,3));
    temp = newobj_move.trialdat(:,:,params.trialid{i});
    newpsth_move(:,:,i) = squeeze(mean(temp,3));
end
newobj.psth = newpsth;
newobj_move.psth = newpsth_move;
clear newpsth newpsth_move

%% calc selectivity and dprime after removing predictions and label cells as sample selective and/or delay selective

cond2use = [2 3]; % right hit - left hit
inp.nTrials = 100;
inp.pval = 0.001;
epoch = 'goCue';


whenSelective = [-1.85 -1.2]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.sample.ts,sel.sample.cluid,shadeix.sample] = calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [-0.8 -0.0]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.delay.ts,sel.delay.cluid,shadeix.delay] = calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

whenSelective = [0.0 0.8]; % seconds relative to epoch, used to find cond a neuron is selective for
[sel.response.ts,sel.response.cluid,shadeix.delay] = calcPrefSelectivity(newobj, params, cond2use, epoch, whenSelective, inp); % preferred - nonpreffered b/w the conds in cond2use

dprime_data = calcdprime(newobj,params,cond2use,inp);
dprime_move = calcdprime(newobj_move,params,cond2use,inp);



%% PLOT DPRIME

close all

baselineedges = [params.eventTimes.bitStart+0.1 params.eventTimes.sample];
baselineix = findTimeIX(obj.time,baselineedges);
bix = baselineix(1):baselineix(2);

cols = getColors;
cmap = flipud(redblue);
alph = 0.2;
dprimethresh = 0.0;
div = 2;

dprime_move_ = dprime_data * 100;

mn = min(min(dprime_move_));
mx = max(max(dprime_move_));

lab = {'sample','delay','response'};
edges = -1:0.1:1;
f = figure; f.Position = [391         227        1199         678];
f.Renderer = 'painters';
% sort by sample epoch dprime
ax = prettifyPlot(subplot(2,3,1));
hold on;
ix = [-1.85 -1.2]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_move_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_move_,2),dprime_move_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / div;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:37;
ylim([notSelective(end)+1 size(dprime_move_,2)])
title('sample')
ylabel('Units')
ax = prettifyPlot(subplot(2,3,4));
hold on;
temp = dprime_move_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)
ylabel('d-prime')
xlabel('Time from go cue (s)')

% sort by delay epoch dprime
ax = prettifyPlot(subplot(2,3,2));
hold on;
ix = [-0.4 0]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_move_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_move_,2),dprime_move_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / div;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:31;
ylim([notSelective(end)+1 size(dprime_move_,2)])
title('delay')
ax = prettifyPlot(subplot(2,3,5));
hold on;
temp = dprime_move_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < -dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)

% sort by response epoch dprime
ax = prettifyPlot(subplot(2,3,3));
hold on;
ix = [0 0.4]; ix = findTimeIX(obj.time,ix); ix = ix(1):ix(2);
[~,sortix] = sort(mean(dprime_move_(ix,:),1),'descend');
imagesc(obj.time,1:size(dprime_move_,2),dprime_move_(:,sortix)');
colormap(cmap); colorbar; 
clims = max(abs(mn),abs(mx)) / div;
clim([-clims clims])
xlim([params.tmin params.tmax])
plotEventTimes(ax,params.eventTimes)
notSelective = 1:69;
ylim([notSelective(end)+1 size(dprime_move_,2)])
title('response')
ax = prettifyPlot(subplot(2,3,6));
hold on;
temp = dprime_move_(:,notSelective(end)+1:end);
tempr = temp(:,nanmedian(temp(ix,:),1) > dprimethresh);
templ = temp(:,nanmedian(temp(ix,:),1) < -dprimethresh);
mu.r = nanmean(tempr,2);
mu.r = mu.r -  mean(mu.r(bix,:),1);
ci.r = getCI(tempr);
mu.l = nanmean(templ,2);
mu.l = mu.l -  mean(mu.l(bix,:),1);
ci.l = getCI(templ);
shadedErrorBar(obj.time,mu.r,ci.r,{'Color',cols.rhit,'LineWidth',2},alph,ax)
shadedErrorBar(obj.time,mu.l,ci.l,{'Color',cols.lhit,'LineWidth',2},alph,ax)
xlim([params.tmin params.tmax])
% ylim([-1 1])
plotEventTimes(ax,params.eventTimes)

