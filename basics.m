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

for isess = 1:numel(meta)
    tag.nTag(isess) = numel(obj(isess).tag);
    tag.tagcluid{isess} = [obj(isess).tag(:).cluid]; % where tagged units are in obj.clu
    tag.tagcluid_params{isess} = find(ismember(params(isess).cluid,tag.tagcluid{isess})); % where in params.cluid, trialdat, psth
end


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

c1 = [255, 87, 225]./255;
c2 = [0, 208, 105]./255;
cmap_ = createcolormap(256, c1, c2);

align = mode(obj.bp.ev.(params.alignEvent));

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
    plotEventTimes(ax,params.eventTimes);
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

    % % KINEMATIC heatmaps

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
        colormap(ax,cmap_)
        c = colorbar; c.Label.String = feat;
        plotEventTimes(ax,params.eventTimes)
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

    % % KINEMATIC trial-avgs

    [f,ax] = prettyPlot();
    f.Position = [321   223   922   194];
    t = tiledlayout('flow');

    for ifeat = 1:numel(feats)
        feat = feats{ifeat};
        featix = find(ismember(kin(1).featLeg,feat));
        ax = nexttile;
        ax = prettifyAxis(ax);
        hold on;

        for icond = 1:numel(cond)
            kindat = kin(isess).dat(:,trials{icond},featix);
            % kindat = removeOutliers(kindat,3);
            kindat = squeeze(nanmean(kindat,2));

            plot(ax,obj(1).time,kindat,'Color',col{icond},'LineWidth',1.5)
        end

        plotEventTimes(ax,params.eventTimes)
        ylabel(feats{ifeat},'Interpreter','none')
        ax.FontSize = 11;
        xlim(xlims)
    end
    xlabel(t,'Time from go cue (s)')

    if sav
        pth = fullfile(utilspth, 'figs', meta(isess).anm, meta(isess).date);
        fn = ['kinematics-tavg'];
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
        plotEventTimes(axRaster,params.eventTimes);
        plotEventTimes(axPSTH,params.eventTimes);
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


%% ISI

% [nISI, ISIcoord, mx] = calcISI(obj,1,-0.005:0.0001:0.005);
% figure;
% plot(ISIcoord* 1000,nISI);

% [f,ax] = prettyPlot();
% for i = 1:numel(obj.clu{1})
%     cla(ax)
%     coord = ISIcoord(:,i);
%     nISI_ = nISI(:,i);
%     plot(coord*1000,nISI_/sum(nISI_),'linewidth',2)
%     ax = prettifyAxis(ax);
%     ylim([0,1])
%     pause
% end

%% correlation b/w single trial fr and jaw disp
close all
feats = {'jaw_ydisp_view1'};
featix = find(ismember(kin(1).featLeg,feats));

ct = 0;
for isess = 1:numel(meta)
    for i = 1:tag.nTag(isess)
        ct = ct + 1;
        clear cc
        trialdat = zscore(squeeze(obj(isess).trialdat(:,tag.tagcluid_params{isess}(i),:)));
        kdat = zscore(squeeze(kin(isess).dat(:,:,featix)));

        for j = 1:size(trialdat,2)
            corrmat = corrcoef(trialdat(:,j),kdat(:,j));
            cc(j) = corrmat(1,2);
        end
        mu(ct) = nanmean(cc);
        
    end
end
f = figure;
f.Position = [682   681   143   229];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
boxplot(mu)
scatter(1,mu,20,'filled','k')
xlim([0.8 1.2])
ylim([min(mu) - 0.1, max(mu) + 0.2])
ax.XTick = [];
ylabel('Correlation to jaw movement')
%% single trial firing rate and jaw disp heatmaps

close all


unit = 1231; % JPV13 2023-10-03
unitix = find(params.cluid==unit);

trialdat = squeeze(obj.trialdat(:,unitix,:));
% feats = {'motion_energy','jaw_ydisp_view1','tongue_angle'};
feats = {'jaw_ydisp_view1'};
featix = find(ismember(kin(1).featLeg,feats));

f = figure;
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
ax2 = prettifyAxis(subplot(1,2,2));
hold on;

fr_y = 1;
k_y = 1;

nTrials = 5;

xlims = [-2 2];

% cols = linspecer(nTrials,'qualitative');
% cols = linspecer;

c1 = [157,78,221]./255;
c2 = [66, 245, 78]./255;
cmap_ = createcolormap(256, c1, c2);
caxis_ = linspace(-0.3,0.7,256);

trials2plot = sort(randsample(obj.bp.Ntrials,nTrials,false));
for i = 1:numel(trials2plot)
    % fr = normalize(trialdat(:,trials2plot(i)),'range',[0 1]);
    fr = normalize(trialdat(:,trials2plot(i)),'zscore');
    fr = fr + (fr_y * (i-1) + max(fr));
    k = normalize(squeeze(kin.dat(:,trials2plot(i),featix)),'range',[0 1]);
    k = k + (k_y * (i-1));

    cc = corrcoef(fr,k);
    cc_(i) = cc(1,2);
    try
        ccix = findTimeIX(caxis_,cc_(i));
    catch
        if cc_(i) < 0
            ccix = 1;
        else
            cc_(i) = 256;
        end
    end

    plot(ax1,obj(1).time,fr,'linewidth',1.5,'color',cmap_(ccix,:))
    plot(ax2,obj(1).time,k,'linewidth',1.5,'color',cmap_(ccix,:))
end
xlim(ax1,xlims)
xlim(ax2,xlims)
plotEventTimes(ax1,evtimes)
plotEventTimes(ax2,evtimes)
xlabel(ax1,'Time from go cue (s)')
ylabel(ax1,'Trials')
title(ax1,'Norm. firing rate')
title(ax2,'Norm. jaw displacement')

hf = figure('Units','normalized','Position',[0.0583    0.5981    0.2203    0.0815]); 
axc = gca;
colormap(axc,cmap_)
hCB = colorbar('north');
set(axc,'Visible',false)
hCB.Position = [0.15 0.3 0.74 0.4];
hf.Position(4) = 0.1000;
hCB.Ticks = linspace(caxis_(1),caxis_(2),5);


%% for a given unit, plot single trials fr and movement
close all


unit = 1231; % JPV13 2023-10-03
unitix = find(params.cluid==unit);

trialdat = squeeze(obj.trialdat(:,unitix,:));
% feats = {'motion_energy','jaw_ydisp_view1','tongue_angle'};
feats = {'jaw_ydisp_view1'};
featix = find(ismember(kin(1).featLeg,feats));

f = figure;
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
ax2 = prettifyAxis(subplot(1,2,2));
hold on;

fr_y = 1;
k_y = 1;

nTrials = 5;

xlims = [-2 2];

% cols = linspecer(nTrials,'qualitative');
% cols = linspecer;

c1 = [157,78,221]./255;
c2 = [66, 245, 78]./255;
cmap_ = createcolormap(256, c1, c2);
caxis_ = linspace(-0.3,0.7,256);

trials2plot = sort(randsample(obj.bp.Ntrials,nTrials,false));
for i = 1:numel(trials2plot)
    % fr = normalize(trialdat(:,trials2plot(i)),'range',[0 1]);
    fr = normalize(trialdat(:,trials2plot(i)),'zscore');
    fr = fr + (fr_y * (i-1) + max(fr));
    k = normalize(squeeze(kin.dat(:,trials2plot(i),featix)),'range',[0 1]);
    k = k + (k_y * (i-1));

    cc = corrcoef(fr,k);
    cc_(i) = cc(1,2);
    try
        ccix = findTimeIX(caxis_,cc_(i));
    catch
        if cc_(i) < 0
            ccix = 1;
        else
            cc_(i) = 256;
        end
    end

    plot(ax1,obj(1).time,fr,'linewidth',1.5,'color',cmap_(ccix,:))
    plot(ax2,obj(1).time,k,'linewidth',1.5,'color',cmap_(ccix,:))
end


%% mike grant plot (heatmaps)
close all

xlims = [-2 2];

unit = 1231; % JPV13 2023-10-03
unitix = find(params.cluid==unit);

c1 = [255, 87, 225]./255;
c2 = [0, 208, 105]./255;
cmap_ = createcolormap(256, c1, c2);

trialdat = squeeze(obj.trialdat(:,unitix,:));
% trialdat = removeOutliers(trialdat,3);
% feats = {'motion_energy','jaw_ydisp_view1','tongue_angle'};
feats = {'jaw_ydisp_view1'};
featix = find(ismember(kin(1).featLeg,feats));
kdat = kin.dat(:,:,featix);
kdat = kdat - min(min(kdat));
% kdat = removeOutliers(kdat,3);

f = figure;
f.Position = [631   521   281   319];
f.Renderer = 'painters';
ax = prettifyAxis(gca); hold on;
imagesc(obj(1).time,1:size(trialdat,2),mySmooth(trialdat,11,'reflect')')
colormap(ax,'parula')
colorbar
clim([0 100])
xlabel('Time from go cue (s)')
ylabel('Trials')
ax.FontSize = 11;
plotEventTimes(ax,evtimes,'k')
ylim([1 size(trialdat,2)+1])
xlim(xlims)

f = figure;
f.Position = [631   521   281   319];
f.Renderer = 'painters';
ax = prettifyAxis(gca); hold on;
imagesc(obj(1).time,1:size(trialdat,2),kdat')
colormap(ax,'turbo')
colorbar
% clim([0 60])
xlabel('Time from go cue (s)')
ylabel('Trials')
ax.FontSize = 11;
plotEventTimes(ax,evtimes,'k')
ylim([1 size(trialdat,2)+1])
xlim(xlims)


%% selectivity

close all

cond2use = [2 3]; % right hit - left hit
epochs = {'goCue','goCue'};
whenSelective = {[-0.5 -0.01], [0.01 0.5]}; % seconds relative to epoch, used to find cond a neuron is selective for
title_string = {'delay','response'};

presample_ix = [evtimes.sample-0.4 evtimes.sample];
presample_ix = findTimeIX(obj(1).time,presample_ix);
presample_ix = presample_ix(1):presample_ix(2);
for i = 1:numel(whenSelective)

    tag.sel = cell(numel(meta),1);
    for isess = 1:numel(meta)
        [sel,shadeix] = calcPrefSelectivity(obj(isess), params(isess), cond2use, epochs{i}, whenSelective{i}); % preferred - nonpreffered b/w the conds in cond2use
        tag.sel{isess} = sel(:,tag.tagcluid_params{isess});
        nontag.sel{isess} = sel(:,~ismember(1:size(sel,2),tag.tagcluid_params{isess}));
    end
    tag.sel_ = cell2mat(tag.sel');
    nontag.sel_ = cell2mat(nontag.sel);

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
    patchx = [obj(1).time(shadeix(1)),obj(1).time(shadeix(end)),obj(1).time(shadeix(end)),obj(1).time(shadeix(1))];
    patchy = [ylims(1) ylims(1) ylims(2) ylims(2)];
    pch = patch(patchx,patchy,[0.3 0.3 0.3]);
    pch.EdgeColor = 'none';
    pch.FaceColor = [0.3 0.3 0.3];
    pch.FaceAlpha = 0.1;
    plot(xlims,[0 0],'k--')
    xlim([-2 2])
    title(title_string{i})
    xlabel('Time from go cue (s)','fontsize',11)
    ylabel('Selectivity | tag (spks/s)','fontsize',11)


    mu = mean(nontag.sel_,2);
    sem = std(nontag.sel_,[],2) ./ sqrt(size(nontag.sel_,2));
    presamp_mu = mean(mean(tag.sel_(presample_ix,:)));
    mu = mu - presamp_mu;
    [f,ax] = prettyPlot();
    f.Position = [680   647   381   231];
    ax = prettifyAxis(ax,0.1);
    hold on;
    shadedErrorBar(obj(1).time,mu,sem,{'color',[60, 61, 60]/255,'linewidth',2},0.2,ax)
    plotEventTimes(ax,evtimes)
    ylims = ax.YLim;
    xlims = ax.XLim;
    patchx = [obj(1).time(shadeix(1)),obj(1).time(shadeix(end)),obj(1).time(shadeix(end)),obj(1).time(shadeix(1))];
    patchy = [ylims(1) ylims(1) ylims(2) ylims(2)];
    pch = patch(patchx,patchy,[0.3 0.3 0.3]);
    pch.EdgeColor = 'none';
    pch.FaceColor = [0.3 0.3 0.3];
    pch.FaceAlpha = 0.1;
    plot(xlims,[0 0],'k--')
    xlim([-2 2])
    title(title_string{i},'fontsize',9)
    xlabel('Time from go cue (s)','fontsize',11)
    ylabel('Selectivity | non-tag (spks/s)','fontsize',11)

end






%% waveforms

close all

[f,ax] = prettyPlot();
ax.Visible = 'off';
f.Color = [0,0,0];
cols = linspecer(sum(tag.nTag));
ct = 1;
for i = 1:numel(meta)
    n = tag.nTag(i);
    for j = 1:n
        if ismember(ct,[1 3 7 8])
            ct = ct + 1;
            continue
        end
        ax = nexttile;
        thisclu = tag.tagcluid{i}(j);
        wv = obj(i).clu{1}(thisclu).spkWavs;
        if size(wv,1) > 50
            wv = wv(43:43*2,:);
        end
        plot(wv,'color',cols(ct,:))
        ct = ct + 1;
        title(num2str(ct-1))
        ax.Color = [0 0 0];
        ax.Visible = 'off';
    end
end

set(f,'InvertHardcopy','off')


%% PCA

cond2use = [2 3];

tag.psth = cell(numel(meta),1);
for isess = 1:numel(meta)
    thisclu = tag.tagcluid_params{isess};
    tag.psth{isess} = obj(isess).psth(:,thisclu,cond2use);
    tag.trialdat{isess} = permute(obj(isess).trialdat(:,thisclu,:),[1 3 2]);
end
tag.psth_ = cell2mat(tag.psth');

dat = permute(tag.psth_,[1 3 2]); % (time,cond,neuron)
sz = size(dat);
temp = reshape(dat,sz(1)*sz(2),sz(3));

temp = zscore(temp);

ncomp = 10;
[pc,~,~,~,ve] = pca(temp,'NumComponents',ncomp); % (neuron,dim)
proj = tensorprod(dat,pc,3,1); % (time,cond,dim)


close all
[f,ax] = prettyPlot();
t = tiledlayout('flow');
cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;
for idim = 1:ncomp
    ax = nexttile;
    ax = prettifyAxis(ax);
    hold on;
    for icond = 1:numel(cond2use)
        plot(obj(1).time,squeeze(proj(:,icond,idim)),'color',col{icond},'linewidth',1.3)
    end
    title(['PC' num2str(idim) ' | VE=' num2str(round(ve(idim),1))])
    plotEventTimes(ax,evtimes)
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Proj (a.u)')

[f,ax] = prettyPlot();
ax = prettifyAxis(ax,0.1);
hold on;
plot(1:numel(ve),ve,'k','LineWidth',2)
scatter(1:numel(ve),ve,30,'filled','MarkerEdgeColor','flat')
ylabel('%VE')
xlabel('PC')

