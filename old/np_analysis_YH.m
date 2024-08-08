clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'nullspace')));

clc

% % TODO
% - fix getStructVarNames to retrive at least trialTypes and Nlicks
% - add stuff to deleteTaggingTrials
% - time warp


%% PARAMETERS

params = defaultParams();

% set conditions to calculate PSTHs for (and get trial numbers for)
params.condition = {};
params.condition(1)         = {'(hit|miss|no)'};
params.condition(end+1)     = {'hit&Nlicks==2&ev.goCue==min(ev.goCue)'}; %ev.goCue==min(ev.goCue) === not early
params.condition(end+1)     = {'hit&Nlicks==5&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==2&trialTypes==1&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==2&trialTypes==2&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==2&trialTypes==3&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==5&trialTypes==1&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==5&trialTypes==2&ev.goCue==min(ev.goCue)'};
params.condition(end+1)     = {'hit&Nlicks==5&trialTypes==3&ev.goCue==min(ev.goCue)'};

params.condLabel{1}     = 'all';
params.condLabel{end+1} = 'hit_2licks';
params.condLabel{end+1} = 'hit_5licks';
params.condLabel{end+1} = 'hit_2licks_right';
params.condLabel{end+1} = 'hit_2licks_center';
params.condLabel{end+1} = 'hit_2licks_left';
params.condLabel{end+1} = 'hit_5licks_right';
params.condLabel{end+1} = 'hit_5licks_center';
params.condLabel{end+1} = 'hit_5licks_left';

% time from align event to grab data for
params.tmin = -2.32;
params.tmax = 4;

params.qm.perform = 0; % Yujin already curated his data

params.events = {'bitStart','goCue'};

params.smooth = 31;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadYH11(meta,datapth);
% meta = meta(1); % 2023-12-06 L_ALM (probe1)
% meta = meta(2); % 2023-12-06 R_ALM (probe2)

meta = loadYH9(meta,datapth); % 2023-11-22 R_FN


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
me = loadMotionEnergy(obj, meta, params, datapth);
kin = getKinematics(obj, me, params);

% TAGGED UNIT META
tagprobenums = [obj.tag(:).probenum];
usetag = tagprobenums==meta.probe;
thistag = obj.tag(usetag);
tag.nTag = numel(thistag);
tag.cluid = [thistag(:).cluid]; % where tagged units are in obj.clu
tag.cluid = find(ismember(params.cluid,tag.cluid))'; % where in params.cluid, trialdat, psth

% REMOVE TAG OVERLAPS
if params.removeTagOverlap
    for i = 1:numel(obj.tag)
        probenum = obj.tag(i).probenum;
        if  probenum ~= meta.probe
            continue
        end
        overlap = obj.tag(i).tagOverlap > 0.9;
        if sum(overlap) > 0
            % where is it in obj.clu{probenum}
            cluix = find(overlap);
            % where is it not in obj.psth/params.cluid
            paramsix = find(~ismember(params.cluid,cluix));
            obj.psth = obj.psth(:,paramsix,:);
            obj.trialdat = obj.trialdat(:,paramsix,:);
            params.cluid = params.cluid(paramsix);
        end
    end
end


%% MOTION ENERGY
close all

cond2plot = [2,3];
trials = cell2mat(params.trialid(cond2plot)');
for i = 1:numel(cond2plot)-1
    yl(i) = numel(params.trialid{cond2plot(i)});
end

f = figure;
f.Position = [506   399   383   450];
ax = prettifyAxis(gca);
hold on;
toplot = me.data(:,trials);
toplot = removeOutliers(toplot,2);
imagesc(obj.time,1:numel(trials),toplot')
colorbar
plotEventTimes(ax,params.eventTimes,'w')
xlim([params.tmin,params.tmax])
ylim([-0,numel(trials)+0.5])
for i = 1:numel(yl)
    ll = plot(ax.XLim,[yl(i),yl(i)],'--','Color','w','LineWidth',2);
end
title('Motion energy')
xlabel('Time from go cue (s)')
ylabel('R2 and R5 Trials');

f = figure;
f.Position = [506   399   383   450];
ax = prettifyAxis(gca);
hold on;
toplot = me.move(:,trials);
imagesc(obj.time,1:numel(trials),toplot')
plotEventTimes(ax,params.eventTimes,'k')
xlim([params.tmin,params.tmax])
ylim([-0,numel(trials)+0.5])
for i = 1:numel(yl)
    ll = plot(ax.XLim,[yl(i),yl(i)],'--','Color','k','LineWidth',2);
end
colormap(ax,(sky))
title('Movement mask')
xlabel('Time from go cue (s)')
ylabel('R2 and R5 Trials');



%% R2 vs R5 coding direction
clear r2
% edges = [0 1];
% edges = [-2 -1.5];
edges = [2 3];
tix = findTimeIX(obj.time,edges,true);

cond2use = [2,3];

r2.dat = nanmean(obj.trialdat(tix,:,params.trialid{2}),3);
% r2.dat = obj.psth(tix,:,2);
r2.mu = nanmean(r2.dat,1);
r2.std = nanstd(r2.dat,[],1);

r5.dat = nanmean(obj.trialdat(tix,:,params.trialid{3}),3);
% r5.dat = obj.psth(tix,:,3);
r5.mu = nanmean(r5.dat,1);
r5.std = nanstd(r5.dat,[],1);
sd = [r2.std,r5.std];

cd = ((r2.mu-r5.mu))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)


%% NULL SPACE

nNPDims = 3;% nNPDims as last arg if you want to change nDims per subspace [4,6,10,13]

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
% -- input data
% trialdat_zscored = zscore_singleTrialNeuralData(obj);
trialdat_zscored = permute(obj.trialdat, [1 3 2]);

% -- null and potent spaces
cond2use = [1]; % all trials
rez = singleTrial_elsayed_np_YH(trialdat_zscored, obj, me, ...
    params, cond2use, nNPDims); % nNPDims as last arg if you want to change nDims per subspace


disp('DONE')

%% projections and reconstructions

proj.cd = tensorprod(trialdat_zscored,cd,3,2);
proj.null = tensorprod(trialdat_zscored,rez.Qnull,3,1);
proj.potent = tensorprod(trialdat_zscored,rez.Qpotent,3,1);

recon.null = tensorprod(proj.null,rez.Qnull,3,2);
recon.potent = tensorprod(proj.potent,rez.Qpotent,3,2);
recon.cd = tensorprod(proj.cd,cd,3,1);

% project null and potent recons onto cd, and reconstruct
proj.nullcd = tensorprod(recon.null,cd,3,2);
proj.potentcd = tensorprod(recon.potent,cd,3,2);

recon.nullcd = tensorprod(proj.nullcd,cd,3,1);
recon.potentcd = tensorprod(proj.potentcd,cd,3,1);

%% variance explained

X = reshape(trialdat_zscored,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
muX = mean(X,1);
total_var_X = sum(var(X-muX,[],1));

% null
Y = reshape(recon.null,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.null = (total_var_X - residual_var) / total_var_X;

% potent
Y = reshape(recon.potent,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.potent = (total_var_X - residual_var) / total_var_X;

% cd
Y = reshape(recon.cd,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.cd = (total_var_X - residual_var) / total_var_X;

% nullcd
Y = reshape(recon.nullcd,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.nullcd = (total_var_X - residual_var) / total_var_X;

% potentcd
Y = reshape(recon.potentcd,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.potentcd = (total_var_X - residual_var) / total_var_X;

%% plot CD projections

cond2plot = [2,3]; % nlicks==2,nlicks==5
trials = params.trialid(cond2plot);

cols = getColors;
c(1,:) = cols.lick2;
c(2,:) = cols.lick5;

dat = proj.cd;

f = figure;
f.Position = [531   540   420   255];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;

for j = 1:numel(cond2plot)
    % this = nanmean(dat(:,trials{j}),2);
    this = dat(:,trials{j});
    [m,ci,mminci,mplusci] = mean_CI(this,0.975);
    shadedErrorBar(obj.time,m,ci,{'LineWidth',2,'Color',c(j,:)},0.2,ax);
    % plot(obj.time,this,'color',c(j,:),'LineWidth',0.1)
end
plotEventTimes(ax,params.eventTimes,'k',false)
xlim(ax,[params.tmin,params.tmax])
ylabel('Proj (a.u)')
xlabel('Time from go cue (s)')
title('R2-R5 CD')


%% plot NP projections

cond2plot = [2,3]; % nlicks==2,nlicks==5
trials = params.trialid(cond2plot);

cols = getColors;
c(1,:) = cols.lick2;
c(2,:) = cols.lick5;

fns = {'null','potent'};
for ispace = 1:numel(fns)
    space = fns{ispace};
    dat = proj.(space);

    f = figure;
    f.Renderer = 'painters';
    t = tiledlayout('flow');

    for i = 1:size(dat,3)
        thisdim = dat(:,:,i);
        ax = prettifyAxis(nexttile);
        hold on;

        for j = 1:numel(cond2plot)
            this = nanmean(thisdim(:,trials{j}),2);
            plot(obj.time,this,'color',c(j,:),'LineWidth',2)
        end

        plotEventTimes(ax,params.eventTimes,'k',false)
        xlim(ax,[params.tmin,params.tmax])
        title([space ' ' num2str(i)],'fontsize',11)
    end
    ylabel(t,'Proj (a.u)')
    xlabel(t,'Time from go cue (s)')
end


%% plot CD NP projections

cond2plot = [2,3]; % nlicks==2,nlicks==5
trials = params.trialid(cond2plot);

cols = getColors;
c(1,:) = cols.lick2;
c(2,:) = cols.lick5;

fns = {'null','potent'};
for ispace = 1:numel(fns)
    space = fns{ispace};
    dat = proj.([space 'cd']);

    f = figure;
    f.Position = [531   540   420   255];
    f.Renderer = 'painters';
    ax = prettifyAxis(gca);
    hold on;

    for j = 1:numel(cond2plot)
        % this = nanmean(dat(:,trials{j}),2);
        this = dat(:,trials{j});
        [m,ci,mminci,mplusci] = mean_CI(this,0.975);
        shadedErrorBar(obj.time,m,ci,{'LineWidth',2,'Color',c(j,:)},0.2,ax);
        % plot(obj.time,this,'color',c(j,:),'LineWidth',0.1)
    end
    plotEventTimes(ax,params.eventTimes,'k',false)
    xlim(ax,[params.tmin,params.tmax])
    ylabel('Proj (a.u)')
    xlabel('Time from go cue (s)')
    title(['R2-R5 CD ' space])
end


%% plot reconstructions
close all

f = figure;
f.Position = [680   646   767   232];
f.Renderer = 'painters';
ax1 = prettifyAxis(subplot(1,2,1));
hold(ax1,'on');
ax2 = prettifyAxis(subplot(1,2,2));
hold(ax2,'on');

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

plt = 1;
if plt

    cond2plot = [2,3];
    for i = 1:rez.N.dims(3) %rez.N.dims(3)-tag.nTag+1:rez.N.dims(3)
        cla(ax1)
        cla(ax2)
        c(1,:) = cols.rhit;
        c(2,:) = cols.lhit;
        for j = 1:numel(cond2plot)
            trix = params.trialid{cond2plot(j)};
            this = squeeze(mean(trialdat_zscored(:,trix,i),2));
            plot(ax1,obj.time,this,'color',c(j,:),'LineWidth',2)
            plot(ax2,obj.time,this,'color',c(j,:),'LineWidth',2)
        end
        c(1,:) = cols.rhit_aw*1.2;
        c(2,:) = cols.lhit_aw*1.2;
        c(c>1) = 1;
        for j = 1:numel(cond2plot)
            trix = params.trialid{cond2plot(j)};
            thisn = squeeze(mean(recon.null(:,trix,i),2));
            thisp = squeeze(mean(recon.potent(:,trix,i),2));
            plot(ax1,obj.time,thisn,'color',c(j,:),'LineWidth',3)
            plot(ax2,obj.time,thisp,'color',c(j,:),'LineWidth',3)
        end


        plotEventTimes(ax1,params.eventTimes,'k',false)
        plotEventTimes(ax2,params.eventTimes,'k',false)
        xlim(ax1,[params.tmin,params.tmax])
        xlim(ax2,[params.tmin,params.tmax])
        title(ax1,'Null')
        title(ax2,'Potent')

        pause

    end
end

%% ve

fns = {'null','potent'};
for j = 1:numel(fns)
    % single trials neural activity reconstructed from n/p
    rc = recon.(fns{j});

    for k = 1:size(rc,3) % for each cell
        orig = trialdat_zscored(:,:,k); % (time,trials)

        pred = rc(:,:,k); % (time,trials) % ve by recon method
        temp = corrcoef(orig(:),pred(:));
        r2.(fns{j})(k) = temp(1,2).^2;
        % mdl = fitlm(orig(:),pred(:));
        % r2.(fns{j})(k) = mdl.Rsquared.Ordinary;
    end
end

totalr2 = r2.null + r2.potent;

%%

alignment = (r2.null - r2.potent) ./ (r2.null + r2.potent);

cols = getColors;

f = figure;
f.Position = [680   581   410   297];
ax = gca;
ax = prettifyAxis(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;

x = find(binedges>0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;

tagalignment = alignment(end-tag.nTag+1:end);
mx = ax.YLim(2) - 4;

scatter(tagalignment,ones(size(tagalignment))*mx, 50,'k','filled','markeredgecolor','w')
xlabel('subspace alignment')
ylabel('unit count')


%% save alignment data

save(['tagalignment_' meta.anm '_' meta.date '.mat'],'tagalignment')

%% plot null aligned psths

cluix = find(alignment > 0.6 & r2.null>0.2);

cond2plot = [2,3]; % nlicks==2,nlicks==5

cols = getColors;


f = figure;
f.Position = [44          45        1832         948];
f.Renderer = 'painters';
t = tiledlayout('flow');

c(1,:) = cols.lick2;
c(2,:) = cols.lick5;


cond2plot = [2,3];
for i = 1:numel(cluix)
    thisclu = cluix(i);
    thispsth = squeeze(obj.psth(:,thisclu,cond2plot));
    thisve = round(r2.null(thisclu)*100,2);
    ax = prettifyAxis(nexttile);
    hold on;

    for j = 1:numel(cond2plot)
        this = thispsth(:,j);
        plot(obj.time,this,'color',c(j,:),'LineWidth',2)
    end

    plotEventTimes(ax,params.eventTimes,'k',false)
    xlim(ax,[params.tmin,params.tmax])
    title(ax,[num2str(thisve) '%'])
end

title(t,'Null-aligned units')
ylabel(t,'Firing rate')
xlabel(t,'Time from go cue (s)')

%% plot potent aligned psths

cluix = find(alignment < -0.6 & r2.potent>0.2);

cond2plot = [2,3]; % nlicks==2,nlicks==5

cols = getColors;


f = figure;
f.Position = [44          45        1832         948];
f.Renderer = 'painters';
t = tiledlayout('flow');

c(1,:) = cols.lick2;
c(2,:) = cols.lick5;


cond2plot = [2,3];
for i = 1:numel(cluix)
    thisclu = cluix(i);
    thispsth = squeeze(obj.psth(:,thisclu,cond2plot));
    thisve = round(r2.potent(thisclu)*100,2);
    ax = nexttile;
    hold on;

    for j = 1:numel(cond2plot)
        this = thispsth(:,j);
        plot(obj.time,this,'color',c(j,:),'LineWidth',2)
    end

    plotEventTimes(ax,params.eventTimes,'k',false)
    xlim(ax,[params.tmin,params.tmax])
    title(ax,[num2str(thisve) '%'])
end

title(t,'Potent-aligned units')
ylabel(t,'Firing rate')
xlabel(t,'Time from go cue (s)')











