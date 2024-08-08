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

params.behav_only = 0;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
meta = loadJPV13(meta,datapth); % 3 sessions

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

%% OBJVIS
% clc,close all
% addpath(genpath(fullfile(utilspth,'objvis')));
% objvis(meta,params,obj,tag,kin)

%% Predict choice from just jaw for now

% feats = {'jaw_ydisp_view1','jaw_xdisp_view1'};
% feats = {'motion_energy'};
feats = {'jaw_ydisp_view1','jaw_xdisp_view1','motion_energy'};
featmask = ismember(kin.featLeg,feats);
kindat = kin.dat(:,:,featmask); % (time,trials,feats)

% PARAMETERS

p.train = 1; % trial percentage
p.test = 1 - p.train;
p.cond2use = [2 5 3 4]; % rhit(went right),lmiss(went right),lhit(went left,rmiss(went left)

p.epoch = 'goCue';
ix = findTimeIX(obj(1).time,[-0.8 0]); % include delay
p.ix = ix(1):ix(2);

% predict trial type
clear trials

% get trials for each lick trial type
labels = {'lickright','lickleft'};
ct = 1;
for i = 1:2:numel(p.cond2use)
    trials.(labels{ct}) = cat(1,params.trialid{p.cond2use(i)}, params.trialid{p.cond2use(i+1)});
    ct = ct + 1;
end

% % sample equal number of trials across lickright and lickleft
t = balanceAndSplitTrials({trials.lickright,trials.lickleft},[1,2],p.train,p.test);
% t.train{1} = train trials for lickright
% t.train{2} = train trials for lickleft
% similar for t.test

% we can test on all the remaining non-train trials actually (TODO)


% predict choice (resp), get X
clear y X
for i = 1:numel(t.train)
    y.train{i} = i * ones(size(t.train{i}));
    y.test{i} = i * ones(size(t.test{i}));

    X.train{i} = kindat(p.ix,t.train{i},:);
    X.test{i} = kindat(p.ix,t.test{i},:);
end
y.train = cell2mat(y.train'); y.test = cell2mat(y.test');
y.train(y.train==1) = 0; y.test(y.test==1) = 0;
y.train(y.train==2) = 1; y.test(y.test==2) = 1;
y.train = logical(y.train); y.test = logical(y.test);

% get design matrix, kinematics
X.train = squeeze(mean(cat(2,X.train{1},X.train{2}),1));
% [nTime,nTrials,nFeats] = size(X.train);
% X.train = reshape(permute(X.train,[2 1 3]), nTrials,nTime*nFeats);

% X.test = cat(2,X.test{1},X.test{2});
X.test = squeeze(mean(cat(2,X.test{1},X.test{2}),1));
% [nTime,nTrials,nFeats] = size(X.test);
% X.test = reshape(permute(X.test,[2 1 3]), nTrials,nTime*nFeats);

% fit
if isrow(X.train)
    X.train = X.train';
end
mdl = fitglm(X.train,y.train,'Distribution','binomial','Link','logit');

% Make predictions
y_proba = predict(mdl, X.train);

% Convert the predicted probabilities to binary outcomes
y_pred = y_proba > 0.5;

acc = sum(y_pred==y.train)/numel(y.train);
disp(['Accuracy: ' num2str(acc)])


%% separate out into trial types

righttrials = find(y.train==0);
lefttrials = find(y.train==1);


% of these trials, which were correctly predicted
r_vr = find(y_pred(righttrials)==0);
l_vl = find(y_pred(lefttrials)==1);


% of these trials, which were incorrectly predicted
r_vl = find(y_pred(righttrials)==1);
l_vr = find(y_pred(lefttrials)==0);

conds.r_vr = t.train{1}(r_vr);
conds.l_vl = t.train{2}(l_vl);
conds.r_vl = t.train{1}(r_vl);
conds.l_vr = t.train{2}(l_vr);

%% plot jaw for each cond
close all

cols = getColors;
c(1,:) = cols.rhit; % r_vr
c(2,:) = cols.lhit; % l_vl
c(3,:) = cols.rhit_aw; % r_vl
c(4,:) = cols.lhit_aw; % l_vr
ls = {'-','-','-','-'};
% c(3,:) = [0.5,0.5,1];
c(4,:) = [245, 130, 218]./255;

feat = {'jaw_ydisp_view1'};
thisfeatdat = kin.dat(:,:,ismember(kin.featLeg,feat));

f = figure;
f.Position = [468   541   605   333];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
cfns = fieldnames(conds);
for i = 1:numel(cfns)
    thisdat = thisfeatdat(:,conds.(cfns{i})); 
    mu = mean(thisdat,2);
    ci = getCI(thisdat);
    shadedErrorBar(obj.time,mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
end
plotEventTimes(ax,params.eventTimes)
xlim([params.tmin,2])
xlabel('Time from go cue (s)')
ylabel('Jaw displacement (y,pixels)')

% same thing zoomed in
tix = findTimeIX(obj.time,[-1.5 0.05],1);

f = figure;
f.Position = [468   478   405   396];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
cfns = fieldnames(conds);
for i = 1:numel(cfns)
    thisdat = thisfeatdat(tix,conds.(cfns{i}));
    mu = mean(thisdat,2);
    ci = getCI(thisdat);
    shadedErrorBar(obj.time(tix),mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
end
plotEventTimes(ax,params.eventTimes)
xlim([-1.5 0.05])
xlabel('Time from go cue (s)')
ylabel('Jaw displacement (y,pixels)')

%% plot ME
close all

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;
c(3,:) = cols.rhit_aw;
c(4,:) = cols.lhit_aw;
ls = {'-','-','-','-'};
% c(3,:) = [0.5,0.5,1];
c(4,:) = [245, 130, 218]./255;

feat = {'motion_energy'};
thisfeatdat = kin.dat(:,:,ismember(kin.featLeg,feat));

f = figure;
f.Position = [468   541   605   333];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
cfns = fieldnames(conds);
for i = 1:numel(cfns)
    thisdat = thisfeatdat(:,conds.(cfns{i}),1); % just the ydisp
    mu = mean(thisdat,2);
    ci = getCI(thisdat);
    shadedErrorBar(obj.time,mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
end
plotEventTimes(ax,params.eventTimes)
xlim([params.tmin,2])
xlabel('Time from go cue (s)')
ylabel('Motion energy')

% same thing zoomed in
tix = findTimeIX(obj.time,[-1.5 0.05],1);

f = figure;
f.Position = [468   478   405   396];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
cfns = fieldnames(conds);
for i = 1:numel(cfns)
    thisdat = thisfeatdat(tix,conds.(cfns{i}),1); % just the ydisp
    mu = mean(thisdat,2);
    ci = getCI(thisdat);
    shadedErrorBar(obj.time(tix),mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
end
plotEventTimes(ax,params.eventTimes)
xlim([-1.5 0.05])
xlabel('Time from go cue (s)')
ylabel('Motion energy')


%% plot tagged units
close all

for itag = 1:tag.nTag

    thisclu = squeeze(obj.trialdat(:,tag.cluid.obj(itag),:));

    f = figure;
    f.Position = [468   541   605   333];
    f.Renderer = 'painters';
    ax = prettifyAxis(gca);
    hold on;
    cfns = fieldnames(conds);
    for i = 1:numel(cfns)
        thisdat = thisclu(:,conds.(cfns{i}));
        mu = mean(thisdat,2);
        ci = getCI(thisdat);
        shadedErrorBar(obj.time,mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
    end
    plotEventTimes(ax,params.eventTimes)
    xlim([params.tmin,2])
    xlabel('Time from go cue (s)')
    ylabel('spks/sec')
    title(['Tagged Unit ' num2str(tag.cluid.clu(itag))])


    % % same thing zoomed in
    % tix = findTimeIX(obj.time,[-1.5 0.05],1);
    % 
    % f = figure;
    % f.Position = [468   478   405   396];
    % f.Renderer = 'painters';
    % ax = prettifyAxis(gca);
    % hold on;
    % cfns = fieldnames(conds);
    % for i = 1:numel(cfns)
    %     thisdat = thisclu(tix,conds.(cfns{i}));
    %     mu = mean(thisdat,2);
    %     ci = getCI(thisdat);
    %     shadedErrorBar(obj.time(tix),mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
    % end
    % plotEventTimes(ax,params.eventTimes)
    % xlim([-1.5 0.05])
    % xlabel('Time from go cue (s)')
    % ylabel('spks/sex')
    % title(['Tagged Unit ' num2str(tag.cluid.clu(itag))])

    
end


%% plot all but tagged units
close all

f = figure;
f.Position = [468   541   605   333];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
for iunit = 1:(numel(params.cluid)-tag.nTag)
    
    cla(ax)

    thisclu = squeeze(obj.trialdat(:,iunit,:));

   
    cfns = fieldnames(conds);
    for i = 1:numel(cfns)
        thisdat = thisclu(:,conds.(cfns{i}));
        mu = mean(thisdat,2);
        ci = getCI(thisdat);
        shadedErrorBar(obj.time,mu,ci,{'Color',c(i,:),'LineWidth',2,'LineStyle',ls{i}},0.2,ax)
    end
    plotEventTimes(ax,params.eventTimes,'k',0 ...
        )
    xlim([params.tmin,2])
    xlabel('Time from go cue (s)')
    ylabel('spks/sec')
    title(['Unit ' num2str(params.cluid(iunit))])

    pause

    
end



