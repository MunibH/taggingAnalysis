clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% Sections
% 1) Load Tag Data
% 2) Plot all PSTHs in one figure
% 3) Plot each PSTH in separate figures
% 4) Plot PCs by region
% 5)

%% LOAD TAG DATA

fpth = 'C:\Users\munib\Documents\Economo-Lab\data\tagged';
% fn = 'AllTagged_GoCue_20240817.mat';
% fn = 'AllTagged_FirstLick_20240817.mat';
fn = 'AllTagged_LastLick_20240817.mat';

load(fullfile(fpth,fn)); % loads 'tag' and 'params'

%% Plot all PSTHs in one figure
close all

regions = {'ALM','M1'}; % can be contained in tag.region

cond = [2,3]; % see tag.condLabel for what each condition is

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

almcol = [255, 119, 0]./255;
tjm1col = [153, 0, 255]./255;

plotError = false;
plotUnitTitle = false;

lw = 1.5;
alph = 0.2;
patchalpha = 0.1;
titlefs = 9;

sm = 11;

xl = [-2.1,params.tmax];

f = figure;
f.Position = [1          41        1920         963];
f.Renderer = 'painters';
t = tiledlayout('flow');

for iregion = 1:numel(regions)
    current_region = regions{iregion};
    for isess = 1:numel(tag) % loop through sessions
        thistag = tag(isess);
        for iprobe = 1:numel(thistag.psth) % loop through probes
            thisregion = thistag.region{iprobe}{1};
            if ~contains(thisregion,current_region)
                continue
            end
            thisdat = thistag.trialdat{iprobe};
            for iunit = 1:size(thisdat,2)
                % ax = prettifyAxis(nexttile);
                ax = nexttile;
                ax.LineWidth = 1;
                hold on;
                for icond = 1:numel(cond)
                    trix = thistag.trialid{cond(icond)};
                    temp = squeeze(thisdat(:,iunit,trix));
                    temp = mySmooth(temp,sm,'reflect');
                    [m,h,lowerbnd,upperbnd] = mean_CI(temp);
                    if plotError
                        shadedErrorBar(thistag.time,m,h,{'Color',c(icond,:),'LineWidth',lw},alph,ax)
                    else
                        plot(ax,thistag.time,m,'Color',c(icond,:),'LineWidth',lw)
                    end
                end
                xlim(xl)
                plotEventTimes(ax,thistag.eventTimes)
                if plotUnitTitle
                    title([thistag.anm ' ' strrep(thistag.date,'-','') ' ' ...
                        thisregion ' U' num2str(thistag.id.clu{iprobe}(iunit))], ...
                        'Interpreter','none', 'FontWeight','normal','FontSize',titlefs)
                end
                if contains(thisregion,'ALM')
                    pcol = almcol;
                elseif contains(thisregion,'M1')
                    pcol = tjm1col;
                end
                p = patch(ax,[ax.XLim flip(ax.XLim)],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],almcol);
                p.FaceAlpha = patchalpha;
                p.FaceColor = pcol;
                p.EdgeAlpha = 0;

                if isess > 1 || iunit > 1
                    ax.XTickLabel = '';
                end
            end
        end
    end
end

xlabel(t,['Time from ' params.alignEvent ' (s)'],'FontSize',18)
ylabel(t,'Spks/sec','FontSize',18)

%% Plot each PSTH in separate figures
close all

cond = [2,3]; % see tag.condLabel for what each condition is

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

almcol = [255, 119, 0]./255;
tjm1col = [153, 0, 255]./255;


lw = 1.5;
alph = 0.2;
patchalpha = 0.1;
titlefs = 13;

xl = [-2.1,params.tmax];


for isess = 1:numel(tag) % loop through sessions
    thistag = tag(isess);
    for iprobe = 1:numel(thistag.psth) % loop through probes

        thisdat = thistag.trialdat{iprobe};
        for iunit = 1:size(thisdat,2)
            f = figure;
            f.Renderer = 'painters';
            f.Position = [695   519   373   229];
            ax = prettifyAxis(gca);
            ax.LineWidth = 1;
            hold on;
            for icond = 1:numel(cond)
                trix = thistag.trialid{cond(icond)};
                temp = squeeze(thisdat(:,iunit,trix));
                [m,h,lowerbnd,upperbnd] = mean_CI(temp);
                shadedErrorBar(thistag.time,m,h,{'Color',c(icond,:),'LineWidth',lw},alph,ax)
            end
            xlim(xl)
            plotEventTimes(ax,thistag.eventTimes)
            thisregion = thistag.region{iprobe}{iunit};
            title([thistag.anm ' ' strrep(thistag.date,'-','') ' ' ...
                thisregion ' U' num2str(thistag.id.clu{iprobe}(iunit))], ...
                'Interpreter','none', 'FontWeight','normal','FontSize',titlefs)
            if contains(thisregion,'ALM')
                pcol = almcol;
            elseif contains(thisregion,'M1')
                pcol = tjm1col;
            end
            p = patch(ax,[ax.XLim flip(ax.XLim)],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],almcol);
            p.FaceAlpha = patchalpha;
            p.FaceColor = pcol;
            p.EdgeAlpha = 0;

            if isess > 1 || iunit > 1
                ax.XTickLabel = '';
            end
            xlabel(['Time from ' params.alignEvent ' (s)'])
            ylabel('Spks/sec')
            drawnow
        end
    end
end

%% Plot PCs by region
close all

% region = 'ALM';
region = 'tjM1';

cond = [2,3]; % see tag.condLabel for what each condition is

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

almcol = [255, 119, 0]./255;
tjm1col = [153, 0, 255]./255;
if contains(region,'ALM')
    pcol = almcol;
elseif contains(region,'M1')
    pcol = tjm1col;
end


lw = 1.5;
alph = 0.2;
patchalpha = 0.1;
titlefs = 13;

xl = [-2.1,params.tmax];

% get PCs and project PSTHs (single trials?)
allpsth = nan(numel(tag(1).time)*numel(cond),200); % (time*nCond,units)
ii = 1;
for isess = 1:numel(tag) % loop through sessions
    thistag = tag(isess);
    for iprobe = 1:numel(thistag.psth) % loop through probes
        thisregion = thistag.region{iprobe}{1};
        if ~contains(thisregion,region)
            continue
        end

        psth = permute(thistag.psth{iprobe}(:,:,cond),[1 3 2]); % (time,cond,units);
        psth = reshape(psth,size(psth,1)*size(psth,2),size(psth,3));
        nUnits = size(psth,2);
        allpsth(:,ii:ii+nUnits-1) = psth;
        ii = ii + nUnits;
    end
end
allpsth(:,ii:end) = []; % remove unused allocation

% PCA
[pc,~,~,~,ve] = pca(allpsth);
cve = cumsum(ve);

% project data onto top 3
proj = allpsth * pc;
proj = proj(:,1:3);

% reshape
proj = reshape(proj,numel(tag(1).time),numel(cond),3); % (time,cond,pc)

% plot VE and proj

f = figure;
f.Renderer = 'painters';
f.Position = [508         526        1244         224];
t = tiledlayout('flow');

ax = prettifyAxis(nexttile);
hold on;
plot(cve,'Color',pcol)
scatter(1:numel(cve),cve,20,'MarkerFaceColor',pcol,'MarkerEdgeColor','none')
plot(ax.XLim,[cve(3),cve(3)],'k--')
xlabel('Cumulative VE (%)')
ylabel('PC #')
title(['Top 3 PCs VE: ' num2str(round(cve(3),2))],'FontSize',titlefs,'FontWeight','normal')

for ipc = 1:size(proj,3)
    ax = prettifyAxis(nexttile);
    hold on;

    thisdat = proj(:,:,ipc);
    for icond = 1:numel(cond)
        temp = squeeze(thisdat(:,icond));
        plot(tag(1).time,temp,'Color',c(icond,:),'LineWidth',lw)
    end
    xlim(xl)
    plotEventTimes(ax,thistag.eventTimes)
    p = patch(ax,[ax.XLim flip(ax.XLim)],[ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],almcol);
    p.FaceAlpha = patchalpha;
    p.FaceColor = pcol;
    p.EdgeAlpha = 0;
    if ipc ==1
        xlabel(ax,['Time from ' params.alignEvent ' (s)'],'FontSize',14)
        ylabel(ax,'Proj (a.u)','FontSize',14)
    end
    title(['PC' num2str(ipc)] ,'FontWeight','normal','FontSize',titlefs)
end
% title(t,['Top 3 PCs for ' region], 'Interpreter','none', 'FontWeight','normal','FontSize',titlefs)

% 3d plot
rez.proj = permute(proj,[1 3 2]); % (time,pc,cond)
fns = fieldnames(tag(1).eventTimes);
for i = 1:numel(fns)
    rez.eventIX(i) = findTimeIX(tag(1).time,tag(1).eventTimes.(fns{i}));
end
rez.ve = ve(1:3);
plt.cond2plot = [1,2];
plt.dim2plot3 = 1:3;
plt.sm = 51;
plt.smtype = 'reflect';
plt.lw = 2;
plt.ms = 40;
plt.col{1} = cols.rhit;
plt.col{2} = cols.lhit;

plotProj3D(rez,plt)







