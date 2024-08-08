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

meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions

% subset meta

% meta = subsetMetaByParams(meta,params);

%%  Decode PT activity from movements
close all

par.sav = 0; % save figure to .fig and .png
par.xlims = [-2 2];
par.cond = [2,3];
par.condLabel = {'rhit','lhit'};
par.sm = 15;
par.rasterCondPad = 10;

par.plotLickRaster = 0;
par.plotPSTH = 1;
par.plotKin = 0;

cols = getColors;
col1{1} = cols.rhit;
col1{2} = cols.lhit;
col2{1} = cols.rhit_aw;
col2{2} = cols.lhit_aw;
par.alph = 0.2;

c1 = [255, 87, 225]./255;
c2 = [0, 208, 105]./255;
cmap_ = createcolormap(256, c1, c2);

par.feats = {'motion_energy','jaw_ydisp_view1','tongue_angle','tongue_length',...
    'tongue_ydisp_view1','nose_ydisp_view1','top_tongue_ydisp_view2','jaw_ydisp_view2'};

par.pre=10; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params.dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% data sets
par.train = 1; % fraction of trials
par.test = 1 - par.train;

par.epochs = {[-1.85,-1.2],[-0.5 0],[0 0.5],[1.5 2]};
par.epochLabels = {'sample','delay','response','outcome'};

rez = struct(); % store results to use across all sessions here

rezct = 1;
for isess = 1:numel(meta)
    clearvars -except datapth meta params utilspth isess par col1 col2 rez rezct

    thismeta = meta(isess);

    % LOAD DATA
    [obj,params] = loadSessionData(thismeta,params);
    me = loadMotionEnergy(obj, thismeta, params, datapth);
    kin = getKinematics(obj, me, params);

    for iepoch = 1:numel(par.epochs)
        par.epochix{iepoch} = findTimeIX(obj.time,par.epochs{iepoch},1);
    end

    % TAGGED UNIT META
    tag.nTag = numel(obj.tag);
    tag.cluid.clu = [obj.tag(:).cluid]; % where tagged units are in obj.clu
    tag.cluid.obj = find(ismember(params.cluid,tag.cluid.clu))'; % where in params.cluid, trialdat, psth

    align = mode(obj.bp.ev.(params.alignEvent));

    trials = params.trialid(par.cond);
    alltrials = cell2mat(trials');

    ifeats = ismember(kin.featLeg,par.feats);
    kindat = kin.dat(:,:,ifeats);

    [nTime,nTrials,~] = size(kindat);

    %%
    for itag = 1:tag.nTag
        thisclunum = tag.cluid.clu(itag);

        tagdat = permute(obj.trialdat(:,tag.cluid.obj(itag),:),[1,3,2]); % (time,trials);

        yout = decodeNeuralFromKin_v2(obj,thismeta,params,par,trials,alltrials,kindat,tagdat);


        % plot psth and predicted psth
        f = figure;
        f.Position = [680   633   339   245];
        ax = prettifyAxis(gca);
        hold on;
        for icond = 1:numel(trials)
            this = mySmooth(yout.input(:,trials{icond}),31,'reflect');
            mu1 = mean(this,2);
            sd = std(this,[],2)./numel(trials{icond}); % omit error bars for now
            shadedErrorBar(obj.time,mu1,sd,{'LineWidth',2,'Color',col1{icond}},par.alph,ax)
            this = yout.pred(:,trials{icond});
            mu2 = mean(this,2);
            sd = std(this,[],2)./numel(trials{icond});
            shadedErrorBar(obj.time,mu2,sd,{'LineStyle','--','LineWidth',3,'Color',col2{icond}},par.alph,ax)
            cctemp = corrcoef(mu1,mu2);
            cc(icond) = cctemp(1,2).^2;
        end
        rez(rezct).r2 = sum(cc)./2;
        plotEventTimes(ax,params.eventTimes);
        xlim(ax,par.xlims)
        ylabel(ax,'Firing rate (spks/s)')
        xlabel(ax,'Time from go cue (s)')
        region = thismeta.region;
        title(ax,[thismeta.anm '' thismeta.date '' region ', Unit' num2str(thisclunum)],'fontsize',10,'FontWeight','normal','Interpreter','none');

        if par.sav
            pth = fullfile(utilspth, 'figs', 'DecodeTaggedUnitsFromMovement');
            fn = [thismeta.anm '_' thismeta.date '_Unit_' num2str(thisclunum) '_' region '_PSTH'];
            mysavefig(f,pth,fn)
        end

        % plot selectivity and predicted selectivity
        f = figure;
        f.Position = [680   633   339   245];
        ax = prettifyAxis(gca);
        hold on;
        sel = mySmooth(mean(yout.input(:,trials{1}),2) - mean(yout.input(:,trials{2}),2), 31, 'reflect');
        this1 = mean(yout.input(:,trials{1}),2) - mean(yout.pred(:,trials{1}),2);
        this2 = mean(yout.input(:,trials{2}),2) - mean(yout.pred(:,trials{2}),2);
        sel_ = mySmooth(this1 - this2, 31, 'reflect');
        plot(obj.time,sel,'LineWidth',2,'Color',[245, 0, 224]./255);
        plot(obj.time,sel_,'LineWidth',2,'Color','k');
        ll = line(ax.XLim,[0 0]);
        ll.LineStyle = '--';
        ll.Color = 'k';
        plotEventTimes(ax,params.eventTimes);
        xlim(ax,par.xlims)
        ylabel(ax,'Sel, $\Delta$ Sel (spks/s)','Interpreter','latex')
        xlabel(ax,'Time from go cue (s)')
        region = thismeta.region;
        title(ax,[thismeta.anm '' thismeta.date '' region ', Unit' num2str(thisclunum)],'fontsize',10,'FontWeight','normal','Interpreter','none');

        if par.sav
            pth = fullfile(utilspth, 'figs', 'DecodeTaggedUnitsFromMovement');
            fn = [thismeta.anm '_' thismeta.date '_Unit_' num2str(thisclunum) '_' region '_Selectivity'];
            mysavefig(f,pth,fn)
        end

        % calc avg selectivity per epoch (pred and actual)
        sel = zscore(mean(yout.input(:,trials{1}),2) - mean(yout.pred(:,trials{2}),2));
        this1 = mean(yout.input(:,trials{1}),2) - mean(yout.pred(:,trials{1}),2);
        this2 = mean(yout.input(:,trials{2}),2) - mean(yout.pred(:,trials{2}),2);
        sel_ = zscore(this1 - this2);
        for iepoch = 1:numel(par.epochs)
            rez(rezct).thissel.(par.epochLabels{iepoch}) = mean(sel(par.epochix{iepoch}));
            rez(rezct).thissel_.(par.epochLabels{iepoch}) = mean(sel_(par.epochix{iepoch}));
        end


        drawnow;
        rezct = rezct + 1;
    end
    %%
end

%% r-squared for tagged units

allr2 = [rez(:).r2];

f = figure;
f.Position = [550   490   303   315];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;

xs = 1;
b = boxchart(allr2);
b.WhiskerLineColor = 'w';
for i = numel(xs)
    this = allr2;
    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(this), 0.25);
    scatter(xx, this, 20,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
end
% xlim([0 2])
ylabel('R-squared')
ylim([0,1])



%%

clear allsel allsel_

lincols = linspecer(numel(par.epochs));

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
for j = 1:numel(par.epochs)
    ax = prettifyAxis(nexttile);
    hold on;

    clear allsel allsel_
    for i = 1:numel(rez)
        allsel(i) = rez(i).thissel.(par.epochLabels{j});
        allsel_(i) = rez(i).thissel_.(par.epochLabels{j});
    end
    scatter(allsel,allsel_,'filled','o','MarkerEdgeColor','none','MarkerFaceColor',lincols(j,:))
    title(par.epochLabels{j})
    axis(ax,'square')
    lims.x = ax.XLim;
    lims.y = ax.YLim;
    mn = min(lims.x(1),lims.y(1));
    mx = min(lims.x(2),lims.y(2));
    plot([-10 10],[-10 10],'k--')
    ax.XLim = [mn mx];
    ax.YLim = [mn mx];
end



%%











