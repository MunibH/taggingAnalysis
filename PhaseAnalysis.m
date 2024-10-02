clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - label licks with which lick number in bout
% - find signifcantly modulated cells during licking/one lick...
%   --- maybe better, find preferred phase for each unit, then sort by that
%       in heatmap

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

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
% meta = loadJPV13(meta,datapth); % 3 sessions
meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)


%% LOAD DATA

thismeta = meta(1);

[obj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(obj,sesspar,thismeta);

me = loadMotionEnergy(obj, thismeta, sesspar, datapth);
kin = getKinematics(obj, me, sesspar);

%% old functions
% assignPhases



%% frequency analysis

[freq,Y] = DoJawFFT(obj,400);

%% get lick pos, vel, phase

p.vidFs = 400;

p.presampMs = 125; % number ms to extract prior to each lick peak
p.postsampMs = 125; % number ms to extract post to each lick peak
p.presamp = round(p.presampMs/1000 * p.vidFs);
p.postsamp = round(p.presampMs/1000 * p.vidFs);
p.peri_center = 15; % window to look around to correct lick peak

% Parameters
p.filt.Fs = p.vidFs;  % Sampling frequency (in Hz, replace with your actual sampling rate)
p.filt.Fc1 = 0.1;    % Lower cutoff frequency (in Hz)
p.filt.Fc2 = 15;   % Upper cutoff frequency (in Hz)

p.fillmiss_window = 10;

[lick] = GetLickPhase(obj,p);
spk = GetUnitSpkPhase(obj,sesspar.cluid,lick);

% figure; plot(lick.peak_pro_phase,1:lick.nLicks,'.') % this is the phase at which
% bpod detected lick (corrected using peak jaw position
% mean(lick.phase) % mean phase, should be ~pi

%% get binned spiking data around licks and preferred phase
clear seq seq_phase pref_phase otest_pval

phasedt = 0.5; % radians

for iprobe = 1:numel(spk)
    for iunit = 1:numel(spk{iprobe})
        % (time from lick start, licknum, unit)
        [seq{iprobe}(:,:,iunit),seq_phase] = getSeqLickTriggered(...
            spk{iprobe}(iunit).phase,spk{iprobe}(iunit).licknum,lick.nLicks,phasedt);

        % Assume `phases` is a 1D array containing the phases (in radians) at which a neuron fired spikes
        pref_phase{iprobe}(iunit) = mod(angle(mean(exp(1i * spk{iprobe}(iunit).phase))), 2*pi);

        % test for circular unimodality
        %  [pval z] = circ_rtest(spk{iprobe}(iunit).phase); % for unimodal data
        [otest_pval{iprobe}(iunit)] = circ_otest(sort(spk{iprobe}(iunit).phase)); % for unimodal and axially biomdal date

    end
end

%% plots - visualize lick pos/vel, and phase

close all

% individual lick trajectories overlaid
cm = flipud(thermal);
cm = cm(1:end,:);
nCols = size(cm,1);
query_points = linspace(1, nCols, lick.nLicks);
cm = interp1(1:nCols, cm, query_points);

f = figure;
f.Renderer = 'painters';
f.Position = [401         418        1006         289];
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
for i = 1:lick.nLicks
    tvec = (1:numel(lick.pos{i}))./p.vidFs * 1000;
    patchline(tvec,normalize(lick.pos{i},'range',[-1 1]),'EdgeColor',cm(i,:))
end
ax2 = prettifyAxis(subplot(1,2,2));
hold on;
for i = 1:lick.nLicks
    tvec = (1:numel(lick.pos{i}))./p.vidFs * 1000;
    patchline(tvec,normalize(lick.vel{i},'range',[-1 1]),'EdgeColor',cm(i,:))
end
xlabel(ax1,'Time (ms)')
ylabel(ax1,'Normalized lick position')
ylabel(ax2,'Normalized lick velocity')
xlim(ax1,[0 250]); xlim(ax2,[0 250])
ylim(ax1,[-1.1 1.1]); ylim(ax2,[-1.1 1.1])


% % lick pos and velocity heatmaps, nans are black
tvec = (1:size(lick.pos_ts,1))./p.vidFs * 1000;

durations = cell2mat(cellfun(@numel,lick.pos,'uni',0));
[~,sortix] = sort(durations);

plotpos = normalize(lick.pos_ts,'range',[-1 1]);
plotpos(isnan(plotpos)) = 10;
plotpos = plotpos(:,sortix);

plotvel = normalize(lick.vel_ts,'range',[-1 1]);
plotvel(isnan(plotvel)) = 10;
plotvel = plotvel(:,sortix);

cmm = brbg;

f = figure;
f.Renderer = 'painters';
f.Position = [516   238   795   566];
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
imagesc(tvec,1:lick.nLicks,plotpos')
c1 = colorbar;
cmap = colormap(cmm);  % Replace 'parula' with your desired colormap
cmap_with_black = [cmap; 0 0 0];
colormap(cmap_with_black);
caxis([-1.1 1]);  % Adjust the lower limit to be below the smallest data value


ax2 = prettifyAxis(subplot(1,2,2));
hold on;
imagesc(tvec,1:lick.nLicks,plotvel')
c2 = colorbar;
cmap = colormap(cmm);  % Replace 'parula' with your desired colormap
cmap_with_black = [cmap; 0 0 0];
colormap(cmap_with_black);
caxis([-1.1 1]);  % Adjust the lower limit to be below the smallest data value

xlabel(ax1,'Time (ms)')
ylabel(ax1,'Lick number')
title(ax1,'Lick position')
title(ax2,'Lick velocity')
xlim(ax1,[0 250]); xlim(ax2,[0 250])
ylim(ax1,[-0.5 lick.nLicks + 0.5]); ylim(ax2,[-0.5 lick.nLicks + 0.5])


% % phase of licks
f = figure;
f.Renderer = 'painters';
f.Position = [516   238   795   566];
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
imagesc(lick.phase_ts_x,1:size(lick.vel_phase,2),normalize(lick.pos_phase(:,sortix),'range',[-1 1])')
c1 = colorbar;
cmap = colormap(cmm);  % Replace 'parula' with your desired colormap
plot(ax1,[pi pi],ax1.YLim,'w--','LineWidth',2)

ax2 = prettifyAxis(subplot(1,2,2));
hold on;
imagesc(lick.phase_ts_x,1:size(lick.vel_phase,2),normalize(lick.vel_phase(:,sortix),'range',[-1 1])')
c2 = colorbar;
cmap = colormap(cmm);  % Replace 'parula' with your desired colormap
plot(ax2,[pi pi],ax2.YLim,'k--','LineWidth',2)

xlabel(ax1,'Phase (rad)')
ylabel(ax1,'Lick number')
title(ax1,'Lick position')
title(ax2,'Lick velocity')
xlim(ax1,[0 2*pi]); xlim(ax2,[0 2*pi])
ylim(ax1,[-0.5 lick.nLicks + 0.5]); ylim(ax2,[-0.5 lick.nLicks + 0.5])


%% plots - num spikes per lick (TODO)

%
f = figure;
f.Position = [879   431   264   251];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
xs = 1:size(nSpikesOnLick,2);
for i = 1:numel(xs)
    this = nSpikesOnLick(:,i);
    xx = simple_violin_scatter(xs(i)*ones(size(this)),this,1000,0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
end
ax.XTick = xs;
xticklabels({'Tag1','Tag2','Tag3','Tag4','Tag5'})
ylabel('# of spikes per lick')
ylim([-1,ax.YLim(2)])

%% plots - lick triggered unit activity

close all

xlc = [235, 52, 174]./255;

confidence = 0.95;
pvalthresh = 1;

qual = sesspar.quality;

f = figure;
f.Position = [636   212   912   752];
f.Renderer = 'painters';
axRaster = prettifyAxis(subplot('Position', [0.1, 0.35, 0.45, 0.6]));
axRaster.XTick = [];
hold on;
axPSTH = prettifyAxis(subplot('Position', [0.1, 0.1, 0.45, 0.22]));
hold on;
axPolar = polaraxes('Position', [0.58, 0.3, 0.4, 0.37]);
axPolar.FontSize = 14;
thetaticks(0:45:360)
axPolar.ThetaAxisUnits = 'radians';
hold on;
for iprobe = 1:numel(spk)
    for iunit = numel(spk{iprobe})-4:numel(spk{iprobe})%1:numel(spk{iprobe})
        cla(axRaster); cla(axPSTH); cla(axPolar)

        if ~(otest_pval{iprobe}(iunit) <= pvalthresh)
            continue
        end

        % raster
        scatter(axRaster,spk{iprobe}(iunit).phase,spk{iprobe}(iunit).licknum,'k.')
        plot(axRaster,[pi pi],axRaster.YLim,'--','Color',xlc,'LineWidth',2)

        % psth
        temp = seq{iprobe}(:,:,iunit);
        % temp = normalize(temp);

        a = 1.0 * temp;

        n = sum(all(~isnan(temp)));
        m = nanmean(a,2);
        se = nanstd(a,[],2) / sqrt(n);
        t_critical = tinv((1 + confidence) / 2, n - 1);
        h = se * t_critical;

        [m,h] = mean_CI(temp);
        shadedErrorBar(seq_phase,m,h,{'Color','k','LineWidth',2},0.3,axPSTH);
        plot(axPSTH,[pi pi],axPSTH.YLim, ...
            '--','Color',xlc,'LineWidth',2)
        scatter(axPSTH,seq_phase,m,'filled','MarkerFaceColor','k','MarkerEdgeColor','flat')
        % plot(axPSTH,[pref_phase{iprobe}(iunit) pref_phase{iprobe}(iunit)],axPSTH.YLim, ...
        %     '--','Color',xlc,'LineWidth',2)

        phist = polarhistogram(axPolar,spk{iprobe}(iunit).phase);
        phist.EdgeColor = 'none';
        phist.FaceColor = 'k';

        ylabel(axRaster,'Lick #')
        title(axRaster,['Probe' num2str(iprobe) ', Unit' num2str(sesspar.cluid{iprobe}(iunit)) ', ' qual{iprobe}{iunit}], ...
            'FontSize', 11)
        xlabel(axPSTH,'Phase (rad)')
        ylabel(axPSTH,'spks/radian')
        xlim(axRaster,[0 2*pi]); xlim(axPSTH,[0 2*pi]);
        ylim(axRaster,[0 lick.nLicks])
        title(axPolar,['p=' num2str(round(otest_pval{iprobe}(iunit),4))],'fontweight','normal','fontsize',12)

        pause

    end
end

%% plots - preferred phase sorted heatmap
close all

cm = brbg;
sm = 1;

pvalthresh = 0.05;

for iprobe = 1:numel(spk)
    f = figure;
    f.Position = [680   239   502   639];
    f.Renderer = 'painters';
    ax = prettifyAxis(gca);
    hold on;

    thisdata = seq{iprobe};
    sz = size(thisdata);
    thisdata_re = reshape(thisdata,sz(1)*sz(2),sz(3));
    basemu = obj.baseline{iprobe}.mu;
    basemu = repmat(basemu',sz(1)*sz(2),1);
    basesd = obj.baseline{iprobe}.sigma;
    basesd = repmat(basesd',sz(1)*sz(2),1);
    thisdata_re = (thisdata_re - basemu) ./ basesd;
    thisdata_zscored_base = reshape(thisdata_re,sz(1),sz(2),sz(3));
    temp = squeeze(nanmean(thisdata_zscored_base,2));
    use = otest_pval{iprobe}<pvalthresh;
    temp = temp(:,use);
    [~, sortOrder] = sort(pref_phase{iprobe}(use));
    toplot = temp(:,sortOrder);
    imagesc(seq_phase,1:size(temp,2),normalize(toplot,'range',[-1 1])')
    colormap(cm)
    cbar = colorbar;
    plot([pi pi],ax.YLim,'k--','LineWidth',2)
    xlabel('Phase (rad)')
    ylabel('Unit #')
    ylim([0.5 size(temp,2)-0.5])
    xlim([0 2*pi])
    title('zscored spks/radian')
end


%% plots - phase tuning for each unit

close all

savfig = 1;
fpth = 'figs/LickPhaseTuning';
fname = [thismeta.anm '_' thismeta.date];

xlc = [235, 52, 174]./255;

confidence = 0.95;

qual = sesspar.quality;

nPerFig = 16;
nClus = sum(cell2mat(cellfun(@numel,sesspar.cluid,'uni',0)));
nFigs = floor(nClus/nPerFig) + double(rem(nClus,nPerFig)>0);

f = figure;
f.Position = [200          84        1161         829];
f.Renderer = 'painters';
t = tiledlayout(4,4);
title(t,strrep(fname,'_',' '))
figct = 1;
plotct = 1;
for iprobe = 1:numel(spk)
    for iunit = 1:numel(spk{iprobe})
        if plotct == nPerFig
            if savfig
                mysavefig(f,fpth,[fname '_' num2str(figct)]);
            end
            figct = figct + 1;
            f = figure;
            f.Position = [200          84        1161         829];
            f.Renderer = 'painters';
            t = tiledlayout(4,4);
            title(t,strrep(fname,'_',' '))
            plotct = 1;
        end
        ax = nexttile();
        ax.Visible = 'off';
        axPolar = polaraxes(t);
        axPolar.FontSize = 9;
        thetaticks(0:45:360)
        axPolar.ThetaAxisUnits = 'radians';
        axPolar.Layout.Tile = plotct; % tile location


        phist = polarhistogram(axPolar,spk{iprobe}(iunit).phase,15,'Normalization','probability');
        % Hide the radial axis (r) and grid
        % axPolar.RAxis.Visible = 'off'; 
        axPolar.ThetaColor = 'k';       % Make the theta (angle) labels black
        % axPolar.RColor = 'none';        % Hide the radial grid lines
        axPolar.ThetaTick = 0:45:360;   % Set the angle ticks (in degrees)

        title(axPolar,[strrep(thismeta.region{iprobe},'_',' ') ', U' num2str(sesspar.cluid{iprobe}(iunit)) ...
            ', p=' num2str(round(otest_pval{iprobe}(iunit),3)) ', ' qual{iprobe}{iunit}],'fontweight','normal','fontsize',10,...
            'Interpreter','none')

        plotct = plotct + 1;
    end
end












