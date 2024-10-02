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

params.smooth = 51;

params.dt = 1/75;


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

%% load motion svd facemap results and ANN neural activity predictions


% load neural activity, motion svd, and neural activity predictions from facemap-ANN
load('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data\MAH24_2024-06-11_FacemapPredictions.mat')
% load('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data\JPV11_2023-06-16_FacemapPredictions.mat')
% pred size = (time*trials,units)
% neuralActivity size = (time*trials,units)
% motionSVD size = (time*trials,SVs)
clearvars -except motionSVD obj sesspar meta thismeta tag params datapth pred ...
    neuralActivity dt nTimeEachTrial epoch_lr epoch_train_loss epoch_ve ve_neuron ...
    me kin


%% put data into trials

nCumTime = cumsum(double(nTimeEachTrial));
for i = 1:obj.bp.Ntrials
    if i == 1
        ix = 1:nCumTime(i);
    else
        ix = nCumTime(i-1)+1:nCumTime(i);
    end
    data.n{i} = neuralActivity(ix,:);
    data.pred{i} = pred(ix,:);
    data.motionSVD{i} = motionSVD(ix,:);

    % get go cue time
    gocue(i) = obj.bp.ev.goCue(i);
end


%% get data aligned to go cue

tmin = params.tmin;
tmax = params.tmax;
cond2plot = [2 3 4 5]; % rhit,lhit   not early

trials = sesspar.trialid(cond2plot);

% find where go cue is in obj.time
gcix = findTimeIX(obj.time,0);
nBeforeGC = (gcix-1);
nAfterGC = numel((gcix+1):numel(obj.time));

nClu = numel(cell2mat(sesspar.cluid'));
nSV = size(data.motionSVD{1},2);

data.plot.n = cell(numel(cond2plot),1);
data.plot.pred = cell(numel(cond2plot),1);
data.plot.svd = cell(numel(cond2plot),1);
for icond = 1:numel(trials)
    thistrials = trials{icond};
    nTrials = numel(thistrials);

    tempN = nan(numel(obj.time),nClu,nTrials);
    tempPred = nan(numel(obj.time),nClu,nTrials);
    tempSVD = nan(numel(obj.time),nSV,nTrials);
    for itrial = 1:nTrials
        thistrial = thistrials(itrial);
        thisgc = gocue(thistrial);

        nTimeThisTrial = size(data.n{thistrial},1);
        thistime = linspace(0,nTimeThisTrial*dt,nTimeThisTrial);
        thisgcix = findTimeIX(thistime,thisgc);
        ix = (thisgcix-nBeforeGC):(thisgcix+nAfterGC);

        nt_ = size(data.n{thistrial},1);
        if ix(end) > nt_
            data.n{thistrial}(end+1:ix(end),:) = nan;
            data.pred{thistrial}(end+1:ix(end),:) = nan;
            data.motionSVD{thistrial}(end+1:ix(end),:) = nan;
        end
        try
            tempN(:,:,itrial) = data.n{thistrial}(ix,:);
        catch
            'a'
        end
        tempPred(:,:,itrial) = data.pred{thistrial}(ix,:);
        tempSVD(:,:,itrial) = data.motionSVD{thistrial}(ix,:);

    end

    data.plot.n{icond} = permute(tempN,[1 3 2]); % (time,trials,neurons)
    data.plot.pred{icond} = permute(tempPred,[1 3 2]);
    data.plot.svd{icond} = permute(tempSVD,[1 3 2]);  % (time,trials,SVs)

end

%% plot training results

close all

f = figure;
f.Renderer = 'painters';
f.Position = [568         593        1221         255];
ax = prettifyAxis(subplot(1,3,1));
hold on;
plot(1:numel(epoch_train_loss),epoch_train_loss,'k','LineWidth',0.1)
scatter(1:numel(epoch_train_loss),epoch_train_loss,20,'filled','MarkerFaceColor','k','MarkerEdgeColor','none')
xlabel('Epoch')
ylabel('Test loss (MSE)')
ax = prettifyAxis(subplot(1,3,2));
hold on;
cc = [121, 199, 125]./255;
plotve = epoch_ve(~(epoch_ve==0));
plot(1:numel(plotve),plotve,'Color',cc,'LineWidth',0.1)
scatter(1:numel(plotve),plotve,30,'filled','MarkerFaceColor',cc,'MarkerEdgeColor','none')
ax.XTickLabel = {'0','75','150'};
ylabel('Variance explained')
ax = prettifyAxis(subplot(1,3,3));
hold on;
cc = [250, 167, 90]./255;
plot(1:numel(epoch_lr),epoch_lr,'Color',cc,'LineWidth',0.1)
scatter(1:numel(epoch_lr),epoch_lr,20,'filled','MarkerFaceColor',cc,'MarkerEdgeColor','none')
ylabel('Learning rate')

%% Plot all PSTHs in figures
close all

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

plotError = true;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

sm = 1;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nCluPerFig = 20;
nFigs = ceil(nClu/nCluPerFig);

for ifig = 1:nFigs

    f = figure;
    f.Position = [1          41        1920         963];
    f.Renderer = 'painters';
    t = tiledlayout('flow');

    clus2plot = (1:nCluPerFig) + ((ifig-1)*nCluPerFig);
    clus2plot(clus2plot>nClu) = [];

    for iunit = 1:numel(clus2plot)
        thisclu = clus2plot(iunit);
        ax = prettifyAxis(nexttile);
        % ax = nexttile;
        % ax.LineWidth = 1;
        hold on;
        for icond = 1:2%numel(cond2plot)
            % plot neural data
            temp = data.plot.n{icond}(:,:,thisclu);
            temp = mySmooth(temp,sm,'reflect');
            [m,h,~,~] = mean_CI(temp);
            if plotError
                shadedErrorBar(obj.time,m,h,{'Color',c(icond,:),'LineWidth',lw},alph,ax)
            else
                plot(ax,obj.time,m,'Color',c(icond,:),'LineWidth',lw)
            end
            % plot prediction
            temp = data.plot.pred{icond}(:,:,thisclu);
            temp = mySmooth(temp,sm,'reflect');
            [m,h,~,~] = mean_CI(temp);
            if plotError
                shadedErrorBar(obj.time,m,h,{'Color',cpred(icond,:),'LineWidth',lwpred},alphpred,ax)
            else
                plot(ax,obj.time,m,'Color',cpred(icond,:),'LineWidth',lwpred)
            end
        end
        xlim(xl)
        plotEventTimes(ax,tag.eventTimes)
        title([qualities{thisclu} ', VE: ' num2str(round(ve_neuron(thisclu)*100,2)) '%'],'fontsize',11.5,'FontWeight','normal')

        if iunit > 1
            ax.XTickLabel = '';
        end
    end
    xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',18)
    ylabel(t,'Spks/sec','FontSize',18)
end

%% Plot tagged unit single trial data
close all

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

lw = 1.5;
lwpred = 1;
alph = 0.1;
alphpred = 0.1;

sm = 15;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nTrials2Plot = 15;
nTrialsCond = cellfun(@(x) size(x,2), data.plot.n,'uni',0);
trialsCond = cellfun(@(x) randsample(x,nTrials2Plot,false), nTrialsCond, 'uni',0);

unitnum = find(ismember(qualities,'tagged'));

f = figure;
f.Renderer = 'painters';
f.Position = [430         391        1299         413];
t = tiledlayout('flow');
for iunit = 1:numel(unitnum)
    ax = prettifyAxis(nexttile);
    hold on;

    thisclu = unitnum(iunit);

    for icond = 1:numel(cond2plot)
        % plot neural data
        temp = data.plot.n{icond}(:,trialsCond{icond},thisclu);
        temp = mySmooth(temp,sm,'reflect');
        plot(ax,obj.time,temp,'Color',c(icond,:),'LineWidth',lw)
        % plot prediction
        temp = data.plot.pred{icond}(:,trialsCond{icond},thisclu);
        plot(ax,obj.time,temp,'Color',cpred(icond,:),'LineWidth',lwpred)
    end

    xlim(xl)
    plotEventTimes(ax,tag.eventTimes)

    if iunit > 1
        ax.XTickLabel = '';
    end


    % title(qualities{thisclu})

end
xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',16)
ylabel(t,'Spks/sec','FontSize',16)

%% Plot tagged unit single trial data as heatmaps
close all

cols = getColors;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

sm = 1;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nTrialsCond = cell2mat(cellfun(@(x) size(x,2), data.plot.n,'uni',0));
nTrialsCond = nTrialsCond(1:end-1);

cmap = inferno;

unitnum = find(ismember(qualities,'tagged'));

for iunit = 1:numel(unitnum)

    thisclu = unitnum(iunit);

    f = figure;
    f.Renderer = 'painters';
    f.Position = [680    47   324   910];
    t = tiledlayout('flow');

    % plot neural data
    ax = nexttile;
    ax.LineWidth = 1;
    ax.FontSize = 12;
    hold on;
    temp = cat(2,data.plot.n{1}(:,:,thisclu),data.plot.n{2}(:,:,thisclu));
    temp = mySmooth(temp,sm,'reflect');
    imagesc(obj.time,1:size(temp,2),temp')
    ylim([0.5,size(temp,2)])
    xlim(xl)
    colormap(cmap)
    colorbar;             % Show the colorbar
    % caxis manual;         % Lock the color axis
    caxis_range = caxis;  % Get the color limits of the first heatmap
    plotEventTimes(ax,tag.eventTimes,'w')
    for i = 1:numel(nTrialsCond)
        plot(ax.XLim,[nTrialsCond(i) nTrialsCond],'w--','LineWidth',2)
    end
    % plot prediction
    ax = nexttile;
    ax.LineWidth = 1;
    ax.FontSize = 12;
    hold on;
    temp = cat(2,data.plot.pred{1}(:,:,thisclu),data.plot.pred{2}(:,:,thisclu));
    temp = mySmooth(temp,sm,'reflect');
    imagesc(obj.time,1:size(temp,2),temp')
    ylim([0.5,size(temp,2)])
    xlim(xl)
    colormap(cmap)
    caxis(caxis_range)
    plotEventTimes(ax,tag.eventTimes,'w')
    for i = 1:numel(nTrialsCond)
        plot(ax.XLim,[nTrialsCond(i) nTrialsCond],'w--','LineWidth',2)
    end

    xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',12)
    ylabel(t,'Trials','FontSize',12)

    % title(qualities{thisclu})

end


%% for tagged units, plot few trials in separate plots

close all
% rng(1)

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

predcol = [237, 109, 235]./255;

sm = 15;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nTrials2Plot = numel(cond2plot);

unitnum = find(ismember(qualities,'tagged'));

for iunit = 1:numel(unitnum)
    f = figure;
    f.Renderer = 'painters';
    f.Position = [430   118   369   686];
    t = tiledlayout('flow');


    thisclu = unitnum(iunit);
    % thisclu = iunit;

    for itrial = 1:nTrials2Plot
        ax = prettifyAxis(subplot(nTrials2Plot,1,itrial));
        hold on;

        % sample random trial from cond2plot{itrial}
        thistrial = randsample(size(data.plot.n{itrial},2),1,false);
        thistrial = 12;
        if itrial==2
            thistrial = 2;
        end
        if itrial==3
            thistrial = 1;
        end
        if itrial==4
            thistrial = 10;
        end
        thistrialid = sesspar.trialid{cond2plot(itrial)}(thistrial);

        traj = obj.traj{1}(thistrialid).ts(:,2,4);
        ft = obj.traj{1}(thistrialid).frameTimes - 0.5;
        % yyaxis right
        % plot(ax,ft,traj,'Color',[0.4,0.4,0.4],'LineWidth',1)

        % plot neural data
        % yyaxis left
        temp = squeeze(data.plot.n{itrial}(:,thistrial,thisclu));
        temp = mySmooth(temp,sm,'reflect');
        plot(ax,obj.time,temp,'Color','k','LineWidth',lw)

        % plot prediction
        temp = squeeze(data.plot.pred{itrial}(:,thistrial,thisclu));
        plot(ax,obj.time,temp,'Color',predcol,'LineWidth',lwpred,'LineStyle','-')

        % if itrial < nTrials2Plot
        %     ax.XTickLabel = '';
        % end
        title(['Trial ' num2str(itrial)],'fontsize',11,'FontWeight','normal')

        xlim(ax,[-2.1 obj.time(end)])
        plotEventTimes(ax,sesspar.eventTimes)


    end

    xlabel('Time (s)','FontSize',12)
    ylabel('Spks/sec','FontSize',12)

    drawnow;
end




%% for tagged units, plot few trials residuals

close all
% rng(1)

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

predcol = [237, 109, 235]./255;

sm = 15;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nTrials2Plot = numel(cond2plot);

unitnum = find(ismember(qualities,'tagged'));

for iunit = 1:numel(unitnum)
    f = figure;
    f.Renderer = 'painters';
    f.Position = [430   118   369   686];
    t = tiledlayout('flow');


    thisclu = unitnum(iunit);
    % thisclu = iunit;

    for itrial = 1:nTrials2Plot
        ax = prettifyAxis(subplot(nTrials2Plot,1,itrial));
        hold on;

        % sample random trial from cond2plot{itrial}
        thistrial = randsample(size(data.plot.n{itrial},2),1,false);
        thistrial = 12;
        if itrial==2
            thistrial = 2;
        end
        if itrial==3
            thistrial = 1;
        end
        if itrial==4
            thistrial = 10;
        end

        thistrialid = sesspar.trialid{cond2plot(itrial)}(thistrial);

        traj = obj.traj{1}(thistrialid).ts(:,2,4);
        ft = obj.traj{1}(thistrialid).frameTimes - 0.5;

        temp = squeeze(data.plot.n{itrial}(:,thistrial,thisclu));
        temp1 = mySmooth(temp,sm,'reflect');
        % plot(ax,obj.time,temp,'Color','k','LineWidth',lw)
        temp2 = squeeze(data.plot.pred{itrial}(:,thistrial,thisclu));
        % plot(ax,obj.time,temp,'Color',predcol,'LineWidth',lwpred,'LineStyle','-')
        resid = temp1 - temp2;
        plot(ax,obj.time,resid,"Color",[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-');

        % if itrial < nTrials2Plot
        %     ax.XTickLabel = '';
        % end
        title(['Trial ' num2str(itrial)],'fontsize',11,'FontWeight','normal')

        xlim(ax,[-2.1 obj.time(end)])
        plotEventTimes(ax,sesspar.eventTimes)


    end

    xlabel('Time (s)','FontSize',12)
    ylabel('Spks/sec','FontSize',12)

    drawnow;
end

%% for tagged units, plot few trials in separate plots along with jaw

close all
% rng(1)

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

predcol = [237, 109, 235]./255;

sm = 15;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

nTrials2Plot = 5;
rconds = randsample([1 2 3 4],nTrials2Plot,true);

unitnum = find(ismember(qualities,'tagged'));

ijaw = find(ismember(kin.featLeg,'motion_energy'));

for iunit = 1:numel(unitnum)
    f = figure;
    f.Renderer = 'painters';
    f.Position = [430   118   369   686];
    t = tiledlayout('flow');


    thisclu = unitnum(iunit);
    % thisclu = iunit;

    for itrial = 1:nTrials2Plot
        ax = prettifyAxis(subplot(nTrials2Plot,1,itrial));
        hold on;

        nTrialsCond = size(data.plot.n{rconds(itrial)},2);
        thistrial = randsample(nTrialsCond,1,false);

        thistrialid = sesspar.trialid{rconds(itrial)+1}(thistrial);

        jaw = mySmooth(squeeze(kin.dat(:,thistrialid,ijaw)),11,'reflect');
        yyaxis right
        cc = [255, 148, 36]./255;
        set(gca, 'YColor', cc./1.2);
        plot(ax,obj.time,jaw,'Color',cc,'LineWidth',1.4)


        % plot neural data
        yyaxis left
        set(gca, 'YColor', [0 0 0]);
        temp = squeeze(data.plot.n{rconds(itrial)}(:,thistrial,thisclu));
        temp = mySmooth(temp,sm,'reflect');
        plot(ax,obj.time,temp,'Color','k','LineWidth',lw)

        % plot prediction
        temp = squeeze(data.plot.pred{rconds(itrial)}(:,thistrial,thisclu));
        plot(ax,obj.time,temp,'Color',predcol,'LineWidth',lwpred,'LineStyle','-')

        % if itrial < nTrials2Plot
        %     ax.XTickLabel = '';
        % end
        title(['Trial ' num2str(itrial)],'fontsize',11,'FontWeight','normal')

        xlim(ax,[-2.1 obj.time(end)])
        plotEventTimes(ax,sesspar.eventTimes)


    end
    yyaxis right
    ylabel('Jaw pos (px)')
    yyaxis left
    xlabel('Time (s)','FontSize',12)
    ylabel('Spks/sec','FontSize',12)

    drawnow;
end



%% plot sample selectivity, delay selectivity, response selectivity
% for tagged units
close all

clear ix
qualities = cat(1,sesspar.quality{:});
unitnum = find(ismember(qualities,'tagged'));

preix = findTimeIX(obj.time,[sesspar.eventTimes.sample-0.28 sesspar.eventTimes.sample],1);
ix.sample = findTimeIX(obj.time,[sesspar.eventTimes.sample sesspar.eventTimes.delay],1);
ix.delay = findTimeIX(obj.time,[sesspar.eventTimes.delay sesspar.eventTimes.goCue],1);
ix.response = findTimeIX(obj.time,[0.01 0.5],1);
epochs = fieldnames(ix);

condMu = cellfun(@(x) squeeze(mean(x,2)), data.plot.n, 'uni', 0);
selectivity = (condMu{1} - condMu{2});
min_vals = min(selectivity, [], 1); % Minimum value for each column
max_vals = max(selectivity, [], 1); % Maximum value for each column
% selectivity = (selectivity - min_vals) ./ (max_vals - min_vals);
% selectivity = selectivity ./ max_vals;
% selectivity = selectivity - mean(selectivity(preix,:),1);

condMu = cellfun(@(x) squeeze(mean(x,2)), data.plot.pred, 'uni', 0);
selectivity_pred = (condMu{1} - condMu{2});
min_vals = min(selectivity_pred, [], 1); % Minimum value for each column
max_vals = max(selectivity_pred, [], 1); % Maximum value for each column
% selectivity_pred = (selectivity_pred - min_vals) ./ (max_vals - min_vals);
% selectivity_pred = selectivity_pred ./ max_vals;
% selectivity_pred = selectivity_pred - mean(selectivity_pred(preix,:),1);


predcol = [237, 109, 235]./255;

f = figure;
f.Position = [596   474   447   268];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
[m,h] = mean_CI(abs(selectivity));
shadedErrorBar(obj.time,m,h,{'Color','k','LineWidth',3},0.4,ax)
[m,h] = mean_CI(abs(selectivity_pred));
shadedErrorBar(obj.time,m,h,{'Color',predcol./1.2,'LineWidth',3},0.4,ax)
xlim([-2.1 obj.time(end)])
plotEventTimes(ax,sesspar.eventTimes);
xlabel('Time from go cue (s)')
ylabel('Selectivity')


% concordance correlation (https://www.alexejgossmann.com/ccc/)
for i = 1:numel(epochs)
    thisepoch = epochs{i};
    thisix = ix.(thisepoch);
    sel1 = selectivity(thisix,:);
    sel2 =  selectivity_pred(thisix,:);
    seldata.(thisepoch) = mean(sel1,1);
    for j = 1:size(sel1,2)
        % temp = corrcoef(sel1(:,j),sel2(:,j));
        % ccc(i,j) = temp(1,2); % (epochs,units)
        ccc(i,j) = concordance_correlation_coefficient(sel1(:,j),sel2(:,j));
    end
end
f = figure;
f.Position = [943   363   292   335];
f.Renderer = "painters";
hold on;
% cmap_sweep(12,purd)
ax = prettifyAxis(gca);
hold on;
xs = 1:3;
for i = 1:numel(xs)
    this = ccc(i,:);
    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(this)./i, 0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
    scatter(xx(unitnum), this(unitnum), 25,'filled', 'markerfacecolor',predcol, 'markeredgecolor','k')
end
ax.XTick = xs;
xticklabels({'sample','delay','response'})
ylabel('CCC')

% plot tagged selectivity and predictions
unitnum = find(ismember(qualities,'tagged'));

f = figure;
f.Position = [343   194   584   701];
f.Renderer = "painters";
t = tiledlayout('flow');
for i = 1:numel(unitnum)
    ax = prettifyAxis(nexttile);
    hold on;
    plot(obj.time,selectivity(:,unitnum(i)),'Color','k','LineWidth',2)
    plot(obj.time,selectivity_pred(:,unitnum(i)),'Color',predcol,'LineWidth',2)
    plotEventTimes(ax,sesspar.eventTimes);
    if i == 1
        ll = legend(ax,{'data','pred'}); ll.Box = 'off';
    end
    title(['Mean CCC: ' num2str(round(mean(ccc(:,unitnum(i))),2))],'fontsize',10)
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Selectivity')
title(t,'Tagged units')
xlim([-2.1 obj.time(end)])


% deltaSelTag = (selectivity(:,unitnum) - selectivity_pred(:,unitnum));
% nTag = numel(unitnum);
% nBoot = 500;
% deltaSelBoot = [];
% for i = 1:nBoot
%     cix = randsample(size(selectivity,2),nTag,true);
%     delta = selectivity(:,cix) - selectivity_pred(:,cix);
%     % delta = mean(delta,2);
%     deltaSelBoot = cat(2,deltaSelBoot,delta);
% end
% f = figure;
% f.Position = [596   474   447   268];
% f.Renderer = "painters";
% cmap_sweep(12,purd)
% ax = prettifyAxis(gca);
% hold on;
% plot(obj.time,mySmooth(deltaSelTag,21,'reflect'),'LineWidth',1.5)
% [m,h,~,~] = mean_CI(deltaSelBoot.^2,0.975,true);
% shadedErrorBar(obj.time,m,h,{'Color','k','LineWidth',2},0.2,ax)
% xlim([-2.1 obj.time(end)])
% plotEventTimes(ax,sesspar.eventTimes)
% xlabel('Time from go cue (s)')
% ylabel('Delta selectivity')

f = figure;
f.Position = [343   658   310   237];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
scatter(abs(seldata.sample),ccc(1,:),8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
scatter(abs(seldata.sample(unitnum)), ccc(1,unitnum), 25,'filled', 'markerfacecolor',predcol, 'markeredgecolor','k')
xlabel('Selectivity in sample epoch')
ylabel('Sample epoch CCC')

f = figure;
f.Position = [343   658   310   237];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
scatter(abs(seldata.delay),ccc(2,:),8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
scatter(abs(seldata.delay(unitnum)), ccc(2,unitnum), 25,'filled', 'markerfacecolor',predcol, 'markeredgecolor','k')
xlabel('Selectivity in delay epoch')
ylabel('Delay epoch CCC')

f = figure;
f.Position = [343   658   310   237];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
scatter(abs(seldata.response),ccc(3,:),8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
scatter(abs(seldata.response(unitnum)), ccc(3,unitnum), 25,'filled', 'markerfacecolor',predcol, 'markeredgecolor','k')
xlabel('Selectivity in response epoch')
ylabel('Response epoch CCC')

%% VE binned

close all

unitnum = find(ismember(qualities,'tagged'));
nTagged = numel(unitnum);

bin_size = 4;
y_true = cat(2,data.plot.n{:}); % (time,trials,units)
y_pred = cat(2,data.plot.pred{:}); % (time,trials,units)

ve_tagged = variance_explained_binned(y_true(:,:,unitnum), y_pred(:,:,unitnum), bin_size);
ve_tagged = mySmooth(ve_tagged,7,'reflect');

nUnits = size(y_true,3);
nBoot = 100;
ve_boot = [];
for i = 1:nBoot
    runits = randsample(nUnits,nTagged,true);
    temp = mean(variance_explained_binned(y_true(:,:,runits), y_pred(:,:,runits), bin_size),2);
    ve_boot = cat(2, ve_boot, temp);
end

tnew = obj.time;

f = figure;
f.Position = [680   489   666   389];
f.Renderer = "painters";
% cm = cmap_sweep(5,rocket);
ax = prettifyAxis(gca);
hold on;
[m,h,~,~] = mean_CI(ve_tagged);
% shadedErrorBar(tnew,m,h,{'Color',cm(2,:),'LineWidth',3},0.4,ax)
plot(tnew,ve_tagged,'LineWidth',2)

hold on;
[m,h,~,~] = mean_CI(ve_boot,0.95,true);
shadedErrorBar(tnew,m,h,{'Color',[0.3 0.3 0.3],'LineWidth',3},0.2,ax)
% plot(tnew,ve_boot,'LineWidth',2)

xlim([-2.1,obj.time(end)])
plotEventTimes(ax,sesspar.eventTimes)

xlabel('Time from go cue (s)')
ylabel('VE')



%% cd late
close all
clear rproj lproj

conds = [1 2]; % rhit,lhit

tix = findTimeIX(obj.time,[sesspar.eventTimes.goCue-0.6 sesspar.eventTimes.goCue],1);

rdata = data.plot.n{1}; % (time,trials,units)
ldata = data.plot.n{2};

rpsth = squeeze(mean(rdata,2));
lpsth = squeeze(mean(ldata,2));

psth = cat(3,rpsth,lpsth);

cd = calcCD(psth,tix);

rproj.data = tensorprod(rdata,cd,3,1);
lproj.data = tensorprod(ldata,cd,3,1);

cols = getColors;

f = figure;
f.Position = [680   611   445   267];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
[m,h] = mean_CI(rproj.data);
shadedErrorBar(obj.time,m,h,{'Color',cols.rhit,'LineWidth',2},0.3,ax);
[m,h] = mean_CI(lproj.data);
shadedErrorBar(obj.time,m,h,{'Color',cols.lhit,'LineWidth',2},0.3,ax);
xlabel('Time from go cue (s)')
ylabel('Proj (au)')
xlim([-2.1 obj.time(end)])
plotEventTimes(ax,sesspar.eventTimes)
title('CDchoice, Data')

rdata = data.plot.pred{1}; % (time,trials,units)
ldata = data.plot.pred{2};

rpsth = squeeze(mean(rdata,2));
lpsth = squeeze(mean(ldata,2));

psth = cat(3,rpsth,lpsth);

cd = calcCD(psth,tix);

rproj.pred = tensorprod(rdata,cd,3,1);
lproj.pred = tensorprod(ldata,cd,3,1);

cols = getColors;

f = figure;
f.Position = [680   611   445   267];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
[m,h] = mean_CI(rproj.pred);
shadedErrorBar(obj.time,m,h,{'Color',cols.rhit,'LineWidth',2},0.3,ax);
[m,h] = mean_CI(lproj.pred);
shadedErrorBar(obj.time,m,h,{'Color',cols.lhit,'LineWidth',2},0.3,ax);
xlabel('Time from go cue (s)')
ylabel('Proj (au)')
xlim([-2.1 obj.time(end)])
plotEventTimes(ax,sesspar.eventTimes)
title('CDchoice, Predictions')


rdata = data.plot.n{1} - data.plot.pred{1}; % (time,trials,units)
ldata = data.plot.n{2} - data.plot.pred{2};

rpsth = squeeze(mean(rdata,2));
lpsth = squeeze(mean(ldata,2));

psth = cat(3,rpsth,lpsth);

cd = calcCD(psth,tix);

rproj.residual = tensorprod(rdata,cd,3,1);
lproj.residual = tensorprod(ldata,cd,3,1);

cols = getColors;

f = figure;
f.Position = [680   611   445   267];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
[m,h] = mean_CI(rproj.residual);
shadedErrorBar(obj.time,m,h,{'Color',cols.rhit,'LineWidth',2},0.3,ax);
[m,h] = mean_CI(lproj.residual);
shadedErrorBar(obj.time,m,h,{'Color',cols.lhit,'LineWidth',2},0.3,ax);
xlabel('Time from go cue (s)')
ylabel('Proj (au)')
xlim([-2.1 obj.time(end)])
plotEventTimes(ax,sesspar.eventTimes)
title('CDchoice, Residuals')


fns = fieldnames(rproj);
for i = 1:numel(fns)
    sel.(fns{i}) = mean(mean(rproj.(fns{i})(tix,:),1),2) - mean(mean(lproj.(fns{i})(tix,:),1),2);
end


percChangeSelDelay = (sel.data-sel.residual)./sel.data;

%% Plot all PSTHs in figures
close all

cols = getColors;
c(1,:) = cols.rhit_aw;
c(2,:) = cols.lhit_aw;
cpred(1,:) = cols.rmiss;
cpred(2,:) = cols.lmiss;

plotError = true;

lw = 2;
lwpred = 2;
alph = 0.1;
alphpred = 0.1;

sm = 1;

xl = [-2.1,params.tmax];

qualities = cat(1,sesspar.quality{:});

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

for iunit = 1:numel(qualities)
    cla(ax1)
    cla(ax2)
    cla(ax3)

    
    rdata = data.plot.n{1}(:,:,iunit);
    rpsth = mean(rdata,2);
    ldata = data.plot.n{2}(:,:,iunit);
    lpsth = mean(ldata,2);
    
    trialdata = cat(2,rdata,ldata);
    
    plot(ax1,obj.time,rpsth,'Color',cols.rhit,'LineWidth',2)
    plot(ax1,obj.time,lpsth,'Color',cols.lhit,'LineWidth',2)
    plotEventTimes(ax1,sesspar.eventTimes)
    xlim(ax1,[-2.1 obj.time(end)])

    imagesc(ax2,obj.time,1:size(trialdata,2),trialdata')
    plotEventTimes(ax2,sesspar.eventTimes)
    xlim(ax2,[-2.1 obj.time(end)])
    ylim(ax2,[-0.1 size(trialdata,2)+0.1])
    

break
end
% xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',18)
% ylabel(t,'Spks/sec','FontSize',18)
