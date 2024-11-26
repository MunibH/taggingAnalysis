clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% calculate CD choice for tagged units and single units, bootstrapped
% calculate CD choice selectivity matrices and transition latency b/w delay
% and response epochs
% see economo 2018, ED fig 7

%% LOAD DATA

datapth = 'C:\Users\munib\Documents\Economo-Lab\data\curated\';

alignEvent = 'goCue'; % 'goCue' 'lastLick' 'firstLick'
SingleAndMulti = true;

ptnfn = findMostRecentFile(fullfile(datapth,'PTNs'), alignEvent);

if SingleAndMulti
    unfn = findMostRecentFile(fullfile(datapth,'UnidentifiedNeurons'), {alignEvent, 'SingleAndMulti'});
else
    unfn = findMostRecentFile(fullfile(datapth,'UnidentifiedNeurons'), {alignEvent, 'SingleUnit'});
end

ptn = load(ptnfn);
un = load(unfn);

% ptn and un should have same number of entries, each corresponds to a session

%%

clearvars -except alignEvent datapth ptnfn unfn utilspth ptn un

%%
clear cond allqual psth

cond = [2 3];
allqual = {};
psth = [];
allregion = [];
for isess = 1:numel(un.sessobj)
    clear temp
    for iprobe = 1:numel(un.sessobj(isess).psth)
        % concat obj.psth and tag.psth for each probe
        temp{iprobe} = cat(2,un.sessobj(isess).psth{iprobe}(:,:,cond),ptn.tag(isess).psth{iprobe}(:,:,cond));
        allregion = cat(1,allregion, cellstr(repmat(un.meta(isess).region{iprobe},numel(sesspar(isess).quality{iprobe}),1)));
    end
    % concat temp across probes
    if numel(temp)>1
        temp2 = cat(2,temp{1},temp{2});
        psth = cat(2,psth,temp2);
        qual = cat(1,un.sesspar(isess).quality{1},un.sesspar(isess).quality{2});
    else
        psth = cat(2,psth,temp{1});
        qual = un.sesspar(isess).quality{1};
    end

    allqual = cat(1,allqual,qual);
end

%%

ibase = 1:20;
temp = permute(psth,[2 1 3]);
baseFR = nanmean(temp(:,ibase,:),[2 3])';

activityRL = cat(1,psth(:,:,1),psth(:,:,2)) - baseFR; % (time x neurons)

%%

rng(42);
clusterRange = 12;
numBootstraps = 1000;
plotTemplates = 1;
[nnmflabels, consistencyScores, H_final] = doNNMF(activityRL', clusterRange, numBootstraps, plotTemplates);


%% plot cluster templates
close all


clusterLabels = unique(nnmflabels);
nClusters = numel(clusterLabels);
nUnitsPerCluster = accumarray(nnmflabels, 1);


% temp = reshape(H_final,size(H_final,1),size(H_final,2)/2,2);
% temp = permute(temp,[2 1 3]); % (time, template, cond)


cols = getColors;
c{1} = cols.rhit;
c{2} = cols.lhit;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');

for j = 1:nClusters
    ax = prettifyAxis(nexttile,'tl',3.5,'fs',13,'lw',2);
    % ax = nexttile;
    hold on
    
    clumask = nnmflabels==clusterLabels(j);
    temp = squeeze(nanmean(psth(:,clumask,:),2));

    for i = 1:2
        plot(un.sessobj(1).time, temp(:,i), 'Color',c{i}, 'LineWidth',3)
    end
    title(['Cluster ' num2str(j) ', ' num2str(nUnitsPerCluster(j)) ' units'],'fontsize',11,'FontWeight','normal')
    plotEventTimes(ax,un.sesspar(1).eventTimes)
end
xlabel(t,'Time from go cue (s)','fontsize',15)
ylabel(t,'Cluster firing rate','fontsize',15)



[s,sortnmf] = sort(nnmflabels);
sel = (psth(:,:,1) - psth(:,:,2)); 
sorted_sel = sel(:,sortnmf,:);

col_min = min(sorted_sel, [], 1); 
col_max = max(sorted_sel, [], 1); 

sorted_sel = (sorted_sel - col_min) ./ (col_max - col_min);

ibase = 1:20;
baseFR = nanmean(sorted_sel(ibase,:),1);
sorted_sel = sorted_sel - baseFR;

f = figure;
f.Renderer = 'painters';
ax = prettifyAxis(gca,'tl',1,'fs',13,'lw',1);
hold on
imagesc(un.sessobj(1).time, 1:size(sorted_sel,2), mySmooth(sorted_sel,21,'reflect')')
xlim([-2.1 2.5])
ylim([-0.05 size(sorted_sel,2)])
colormap(flipud(balanced))
colorbar
plotEventTimes(ax,un.sesspar(1).eventTimes,'k')
total = 0;
for i = 2:nClusters
    total = total + sum(s==(i-1));
    plot(ax.XLim, [total total],'k-','linewidth',1.5)
    % plot(ax.XLim, [total total],'-','linewidth',2,'Color',[255, 235, 107]./255)
end
title('Selectivity','fontsize',11,'fontweight','normal')
xlabel('Time from go cue (s)', 'fontsize',15)
ylabel('Neurons', 'fontsize',15)
ax.YDir = 'reverse';


%% tagging data

clear cond allqual psth

cond = [1];
allqual = {};
psth = [];
for isess = 1:numel(un.sessobj)
    clear temp
    for iprobe = 1:numel(un.sessobj(isess).psth)
        % concat obj.psth and tag.psth for each probe
        temp{iprobe} = cat(2,un.sessobj(isess).psth{iprobe}(:,:,cond),ptn.tag(isess).psth{iprobe}(:,:,cond));
    end
    % concat temp across probes
    if numel(temp)>1
        temp2 = cat(2,temp{1},temp{2});
        psth = cat(2,psth,temp2);
        qual = cat(1,un.sesspar(isess).quality{1},un.sesspar(isess).quality{2});
    else
        psth = cat(2,psth,temp{1});
        qual = un.sesspar(isess).quality{1};
    end

    allqual = cat(1,allqual,qual);
end


close all


clusterLabels = unique(nnmflabels);
nClusters = numel(clusterLabels);
nUnitsPerCluster = accumarray(nnmflabels, 1);


cols = getColors;
c{1} = cols.rhit;
c{2} = cols.lhit;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');

for j = 1:nClusters
    ax = prettifyAxis(nexttile,'tl',3.5,'fs',13,'lw',2);
    hold on
    
    clumask = nnmflabels==clusterLabels(j);
    temp = squeeze(nanmean(psth(:,clumask),2));
    plot(un.sessobj(1).time, normalize(temp), 'Color',cols.un, 'LineWidth',3)

    clumask = nnmflabels==clusterLabels(j) & ismember(allqual,'tagged');
    temp = nanmean(psth(:,clumask),2);
    plot(un.sessobj(1).time, mySmooth(normalize(temp),21,'reflect'), 'Color',cols.ptn, 'LineWidth',3)


    title(['Cluster ' num2str(j) ', ' num2str(nUnitsPerCluster(j)) ' units'],'fontsize',11,'FontWeight','normal')
    plotEventTimes(ax,un.sesspar(1).eventTimes)
end
xlabel(t,'Time from go cue (s)','fontsize',15)
ylabel(t,'Cluster firing rate','fontsize',15)


f = figure;
f.Renderer = 'painters';
f.Position = [680   632   565   246];
ax = prettifyAxis(gca,'tl',2,'fs',13,'lw',2);
hold on
for j = 1:nClusters
    n = sum(nnmflabels==clusterLabels(j));
    ntag = sum(nnmflabels==clusterLabels(j) & ismember(allqual,'tagged'));
    prop(j) = ntag/n * 100;
end
bar(1:numel(prop), prop)
xlabel('Cluster #')
ylabel('PTNs in cluster (%)')
ax.XTick = 1:nClusters;
% ax.XTickLabel = 

f = figure;
f.Renderer = 'painters';
f.Position = [680   632   565   246];
ax = prettifyAxis(gca,'tl',2,'fs',13,'lw',2);
hold on
for j = 1:nClusters
    clumask = nnmflabels==clusterLabels(j);
    
    meanfr(j) = nanmedian(nanmedian(psth(:,clumask)));
end
bar(1:numel(meanfr), meanfr)
xlabel('Cluster #')
ylabel('Median firing rate')
ax.XTick = 1:nClusters;
% ax.XTickLabel = 


%%

save('C:\Users\munib\Documents\Economo-Lab\data\curated\nnmflabels_Selectivity_AllSingleAndMultiUnit_goCue_20241125.mat','nnmflabels')










