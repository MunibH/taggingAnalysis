clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'rastermap')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath('C:\npy-matlab\npy-matlab'))

clc


%% LOAD TAG DATA

fpth = 'C:\Users\munib\Documents\Economo-Lab\data\curated';

fn = fullfile('PTNs','AllTagged_goCue_20241023.mat');
load(fullfile(fpth,fn)); % loads 'tag' and 'params'

fn = fullfile('UnidentifiedNeurons','AllSingleUnit_goCue_20241023.mat');
load(fullfile(fpth,fn)); % loads '' and 'params'

%%
clear cond allqual psth

cond = [2 3];
allqual = {};
psth = [];
for isess = 1:numel(sessobj)
    clear temp
    for iprobe = 1:numel(sessobj(isess).psth)
        % concat obj.psth and tag.psth for each probe
        temp{iprobe} = cat(2,sessobj(isess).psth{iprobe}(:,:,cond),tag(isess).psth{iprobe}(:,:,cond));
    end
    % concat temp across probes
    if numel(temp)>1
        temp2 = cat(2,temp{1},temp{2});
        psth = cat(2,psth,temp2);
        qual = cat(1,sesspar(isess).quality{1},sesspar(isess).quality{2});
    else
        psth = cat(2,psth,temp{1});
        qual = sesspar(isess).quality{1};
    end
    
    allqual = cat(1,allqual,qual);
end

%%

ibase = 1:20;
temp = permute(psth,[2 1 3]);
baseFR = nanmean(temp(:,ibase,:),[2 3])';

activityRL = cat(1,psth(:,:,1),psth(:,:,2)) - baseFR;

%% run tSNE 10 times, pick the run with the lowest loss

% activityRL is probably a (neurons x time*cond) matrix of PSTHs
nRuns = 1;

mappedX_allRun = {};
loss_tSNE_allRun = [];
for i_run = 1:nRuns
    % [mappedX_iRun loss_tSNE_iRun]=tsne(activityRL','Algorithm','exact','Distance','cosine','NumDimensions',2,'NumPCAComponents',30,'Perplexity',50,'Verbose',2);
    [mappedX_iRun loss_tSNE_iRun]=tsne(activityRL','Algorithm','exact','Distance','cosine',...
        'NumDimensions',2,'NumPCAComponents',50,'Perplexity',50,'Verbose',0);
    mappedX_allRun{i_run,1} = mappedX_iRun;
    loss_tSNE_allRun(i_run,1) = loss_tSNE_iRun;
end
[dummy i_best_run] = min(loss_tSNE_allRun);
mappedX = mappedX_allRun{i_best_run};
clear mappedX_allRun  mappedX_allRun  mappedX_iRun  loss_tSNE_allRun  dummy i_best_run

%%

cols = getColors;

iptn = ismember(allqual,'tagged');
iun = ~iptn;

f = figure; 
f.Position = [458   518   395   309];
ax = prettifyAxis(gca);
hold on

scatter(mappedX(iun,1), mappedX(iun,2), 20, cols.un,'filled');
scatter(mappedX(iptn,1), mappedX(iptn,2), 25, cols.ptn,'filled');
xlabel('Dim 1');
ylabel('Dim 2');
title('t-SNE, all sessions');

%%

ix1 = mappedX(:,1) <= 3;
ix2 = mappedX(:,1) > 3;


f = figure; 
f.Position = [700   457   290   330];

ax = prettifyAxis(subplot(2,1,1));
hold on
plot(sessobj(1).time, nanmean(psth(:,ix1,1),2), 'Color', cols.rhit, 'LineWidth', 2);
plot(sessobj(1).time, nanmean(psth(:,ix1,2),2), 'Color', cols.lhit, 'LineWidth', 2);
xlim([-2 3])
plot([-1.85 -1.85],ax.YLim,'k--')
plot([-1.2 -1.2],ax.YLim,'k--')
plot([0 0],ax.YLim,'k--')
title('Dim1 < 3','fontsize',10.5)

ax = prettifyAxis(subplot(2,1,2));
hold on
plot(sessobj(1).time, nanmean(psth(:,ix2,1),2), 'Color', cols.rhit, 'LineWidth', 2);
plot(sessobj(1).time, nanmean(psth(:,ix2,2),2), 'Color', cols.lhit, 'LineWidth', 2);
xlim([-2 3])
plot([-1.85 -1.85],ax.YLim,'k--')
plot([-1.2 -1.2],ax.YLim,'k--')
plot([0 0],ax.YLim,'k--')
xlabel('Time from go cue (s)');
ylabel('spks/sec');
title('Dim1 > 3','fontsize',10.5)
















