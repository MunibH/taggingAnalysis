clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% % TODO
% - selectivity
% - correlations and xcorrelations with movements
% - kinematic decoding like in our paper
% - np analysis
% - - need to automate finding ME threshold



%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'firstLick';


%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
% meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
% meta = loadJPV13(meta,datapth);

% meta = meta(1);

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



for isess = 1:numel(meta)
    tag.nTag(isess) = numel(obj(isess).tag);
    tag.tagcluid{isess} = [obj(isess).tag(:).cluid]; % where tagged units are in obj.clu
    tag.tagcluid_params{isess} = find(ismember(params(isess).cluid,tag.tagcluid{isess})); % where in params.cluid, trialdat, psth
end

% GET EVENT TIMES
evtimes = getEventTimes(obj(1).bp.ev,params(1).events,params(1).alignEvent);
align = mode(obj(1).bp.ev.(params(1).alignEvent));


%% PLOT LICK RASTER and KINEMATICS and RASTER/PSTH

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

xlims = [-2 2];

cond = [2,3];

sav = 0;

sm = 15;
rasterCondPad = 10;

feats = {'motion_energy','jaw_ydisp_view1','tongue_angle'};

close all

for isess = 1:numel(meta)

    trials = params(isess).trialid(cond);
    alltrials = cell2mat(trials');

    % % LICK RASTER

    [f,ax] = prettyPlot();
    f.Position = [680   560   258   318];
    ax = prettifyAxis(ax);
    hold on;

    lick.l = cellfun(@(x) x-align, obj(isess).bp.ev.lickL,'uni',0);
    lick.r = cellfun(@(x) x-align, obj(isess).bp.ev.lickR,'uni',0);


    lasttrial = 0;
    for j = 1:numel(cond)
        clear allLicks
        if j > 1; pad = rasterCondPad; else; pad = 0; end

        % trix = params(isess).trialid{cond(j)};
        trix = trials{j};
        lick.trialtm.l = lick.l(trix);
        lick.trialtm.r = lick.r(trix);
        lick.trial = renum(trix) + lasttrial + pad;

        lasttrial = lick.trial(end);

        % need to make lick.trialtm.l/r and lick.trial same size

        allLicks.ltm = cell2mat(lick.trialtm.l');
        allLicks.rtm = cell2mat(lick.trialtm.r');

        rindex = 1;
        lindex = 1;
        for itrial = 1:numel(trix)
            nLicksL = numel(lick.trialtm.l{itrial});
            nLicksR = numel(lick.trialtm.r{itrial});
            allLicks.lTrial(lindex:lindex + nLicksL) = lick.trial(itrial);
            allLicks.rTrial(rindex:rindex + nLicksR) = lick.trial(itrial);
            rindex = rindex + nLicksR;
            lindex = lindex + nLicksL;
        end
        allLicks.rTrial = allLicks.rTrial(1:end-1);
        allLicks.lTrial = allLicks.lTrial(1:end-1);


        plot(ax,allLicks.rtm,allLicks.rTrial,'.','color',col{1});
        plot(ax,allLicks.ltm,allLicks.lTrial,'.','color',col{2});

    end
    plotEventTimes(ax,evtimes);
    ax.FontSize = 11;
    xlim(ax,xlims)
    title([meta(isess).anm ' ' meta(isess).date],'fontsize',10)
    ylabel(ax,'Trial')
    xlabel(ax,'Time from go cue (s)')

    if sav
        pth = fullfile(utilspth, 'figs', meta(isess).anm, meta(isess).date);
        fn = ['lickRaster'];
        mysavefig(f,pth,fn)
    end

    % % KINEMATICS

    [f,ax] = prettyPlot();
    f.Position = [398         370        1008         378];
    t = tiledlayout('flow');

    for ifeat = 1:numel(feats)
        feat = feats{ifeat};
        featix = find(ismember(kin(1).featLeg,feat));
        kindat = kin(isess).dat(:,alltrials,featix);
        kindat = removeOutliers(kindat,3);

        ax = nexttile;
        ax = prettifyAxis(ax);
        hold on;


        imagesc(ax,obj(1).time,1:size(kindat,2),kindat','interpolation','bilinear')
        % colormap(ax,linspecer)
        c = colorbar; c.Label.String = strrep(feat,'_','-');
        clim([130 190])
        plotEventTimes(ax,evtimes)
        ax.FontSize = 11;
        xlim(xlims)
        yl = cumsum(cell2mat(cellfun(@(x) numel(x),trials,'uni',0)));
        for i = 1:numel(trials)-1
            yline(yl(i),'k--')
        end
        ax = prettifyAxis(ax);
        ylim([1,numel(alltrials)])

    end
    ylabel(t,'Trials')
    xlabel(t,'Time from go cue (s)')

    if sav
        pth = fullfile(utilspth, 'figs', meta(isess).anm, meta(isess).date);
        fn = ['kinematics'];
        mysavefig(f,pth,fn)
    end


    % % RASTER AND PSTH

    % psth
    for i = 1:tag.nTag(isess)
        [f,ax] = prettyPlot();
        f.Position = [680   560   258   318];
        axRaster = subplot(2,1,1);
        axRaster = prettifyAxis(axRaster);
        axPSTH = subplot(2,1,2);
        axPSTH = prettifyAxis(axPSTH);
        hold(axPSTH,'on');
        hold(axRaster,'on');

        tagdat = mySmooth(squeeze(obj(isess).psth(:,tag.tagcluid_params{isess}(i),cond)),sm);

        lasttrial = 0;
        for j = 1:numel(cond)
            plot(axPSTH,obj(isess).time,tagdat(:,j),'color',col{j},'linewidth',2)

            if j > 1; pad = rasterCondPad; else; pad = 0; end

            % trix = params(isess).trialid{cond(j)};
            trix = trials{j};
            thisclunum = params(isess).cluid(tag.tagcluid_params{isess}(i));
            thisclu = obj(isess).clu{1}(thisclunum);
            mask = ismember(thisclu.trial,trix);
            trial = thisclu.trial(mask);
            trialtm = thisclu.trialtm(mask) - obj(isess).bp.ev.(params(isess).alignEvent)(trial);
            trial = renum(trial) + lasttrial + pad;


            lasttrial = trial(end);

            plot(axRaster,trialtm,trial,'.','color',col{j});

        end
        plotEventTimes(axRaster,evtimes);
        plotEventTimes(axPSTH,evtimes);
        axRaster.FontSize = 11;
        axPSTH.FontSize = 11;
        xlim(axRaster,xlims)
        xlim(axPSTH,xlims)
        title(axRaster,['Unit ' num2str(thisclunum)],'fontsize',10)
        ylabel(axRaster,'Trial')
        ylabel(axPSTH,'Firing rate (spks/s)')
        xlabel(axPSTH,'Time from go cue (s)')


        if sav
            pth = fullfile(utilspth, 'figs', meta(isess).anm, meta(isess).date);
            fn = ['TaggedUnit_' num2str(thisclunum) '_RasterPSTH'];
            mysavefig(f,pth,fn)
        end
    end

end

%% selectivity

close all

cond2use = [2 3]; % right hit - left hit
sm = 11;
nBoot = 100;
ct = 1;
for isess = 1:numel(meta)
    sel = calcPrefSelectivity_EachTimePoint(obj(isess), params(isess), cond2use, sm);
    tag.sel{isess} = sel(:,tag.tagcluid_params{isess});

    nSamp = numel(tag.tagcluid_params{isess});
    tempdat = sel(:,~ismember(1:size(sel,2),tag.tagcluid_params{isess})); % nontag data to sample from
    for iboot = 1:nBoot
        sampix = randsample(size(tempdat,2),nSamp,true);
        nontag.sel(:,ct) = mean(tempdat(:,sampix),2);
        ct = ct + 1;
    end

end
tag.sel_ = cell2mat(tag.sel);

presample_ix = [evtimes.bitStart+0.1 evtimes.sample];
presample_ix = findTimeIX(obj(1).time,presample_ix);
presample_ix = presample_ix(1):presample_ix(2);

mu = mean(tag.sel_,2);
sem = std(tag.sel_,[],2) ./ sqrt(sum(tag.nTag));
presamp_mu = mean(mean(tag.sel_(presample_ix,:)));
mu = mu - presamp_mu;
[f,ax] = prettyPlot();
f.Position = [680   647   381   231];
ax = prettifyAxis(ax,0.1);
hold on;
shadedErrorBar(obj(1).time,mu,sem,{'color',[21, 84, 21]/255,'linewidth',2},0.2,ax)
plotEventTimes(ax,evtimes)
ylims = ax.YLim;
xlims = ax.XLim;
plot(xlims,[0 0],'k--')
xlim([-2 2])
xlabel('Time from go cue (s)','fontsize',11)
ylabel('Selectivity | tag (spks/s)','fontsize',11)

mu = mean(nontag.sel,2);
sem = std(nontag.sel,[],2) ./ sqrt(nBoot); % ./ sqrt(size(nontag.sel,2));
presamp_mu = mean(mean(nontag.sel(presample_ix,:)));
mu = mu - presamp_mu;
[f,ax] = prettyPlot();
f.Position = [680   647   381   231];
ax = prettifyAxis(ax,0.1);
hold on;
shadedErrorBar(obj(1).time,mu,sem,{'color',[0.5 0.5 0.5]/255,'linewidth',2},0.2,ax)
plotEventTimes(ax,evtimes)
ylims2 = ax.YLim;
xlims = ax.XLim;
plot(xlims,[0 0],'k--')
xlim([-2 2])
ylim([ylims(1) ylims2(2)])
xlabel('Time from go cue (s)','fontsize',11)
ylabel('Selectivity | nontag (spks/s)','fontsize',11)

%% movement correlations

clear cc movecorr
if isfield(tag,'kincc')
    tag = rmfield(tag,{'kincc','kincc_'});
end

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

xlims = [-2 2];

cond = [2,3];

sav = 0;

sm = 15;

feats = {'motion_energy','jaw_ydisp_view1','tongue_angle', 'tongue_length','nose_ydisp_view1'};
for isess = 1:numel(meta)
    neuraldat = permute(obj(isess).trialdat,[1 3 2]); % (time,trials,neurons)
    for ifeat = 1:numel(feats)
        feat = feats{ifeat};
        featix = find(ismember(kin(1).featLeg,feat));
        kindat = kin(isess).dat(:,:,featix);
        kindat = removeOutliers(kindat,3);
        kindat = fillmissing(kindat,"constant",0);
        for iclu = 1:size(neuraldat,3)
            neuraldat_iclu = neuraldat(:,:,iclu);
            tempcorr = corrcoef(kindat(:),neuraldat_iclu(:));
            cc{isess}(ifeat,iclu) = tempcorr(1,2);
        end
    end
    tag.kincc{isess} = cc{isess}(:,tag.tagcluid_params{isess});
    nontag.kincc{isess} = cc{isess}(:,~ismember(1:size(cc{isess},2),tag.tagcluid_params{isess}));
end

tag.kincc_ = cell2mat(tag.kincc);
nontag.kincc_ = cell2mat(nontag.kincc);



%%

close all
f = figure;
f.Renderer = 'painters';
ax = prettifyAxis(gca,2,15);
hold on;
xs = 1:numel(feats);
cols = linspecer(numel(feats));
for ifeat = 1:numel(feats)

    this = tag.kincc_(ifeat,:);
    mu = median(this);

    xx = simple_violin_scatter(xs(ifeat)*ones(size(this)), this, numel(meta), 0.4); 
    scatter(xx,this,40,'filled','markeredgecolor','none','markerfacecolor',cols(ifeat,:))

end
ax.XTick = xs;
ff = strrep(feats,'_','-');
xticklabels(ff)
ylabel('Correlation')

% xx = 1:numel(feats);
% scatter(xx,tag.kincc_,'filled','MarkerEdgeColor','w','MarkerFaceColor','g')
% scatter(xx+0.3,nontag.kincc_,'filled','MarkerEdgeColor','w','MarkerFaceColor','k')


%%

close all
[f,ax] = prettyPlot();
t = tiledlayout('flow');
% f.Position = [680   647   381   231];
for ifeat = 1:numel(feats)
    ax = nexttile;
    ax = prettifyAxis(ax);
    hold on;
    histogram(tag.kincc_(ifeat,:),15,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','g','Normalization','probability')
    histogram(nontag.kincc_(ifeat,:),'EdgeColor','none','FaceAlpha',0.5,'FaceColor','k','Normalization','probability')
    title(feats{ifeat},'Interpreter','none')
end
xlabel(t,'Corr')
ylabel('Probability')

% xx = 1:numel(feats);
% scatter(xx,tag.kincc_,'filled','MarkerEdgeColor','w','MarkerFaceColor','g')
% scatter(xx+0.3,nontag.kincc_,'filled','MarkerEdgeColor','w','MarkerFaceColor','k')