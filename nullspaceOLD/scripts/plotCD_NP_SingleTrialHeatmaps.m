%%

cm = redblue;
% cm = parula;
close all

cdix = 1; % choice (JEB7 2021-04-29)
cond2use = [6,5];
% cdix = 3; % ramp (JEB6 2021-04-18)
% cond2use = [9 9];

nullsm = 7;
potentsm = 0;

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay'},'goCue');

% only show delay period
delayedges = [-2.39 0];
delayt = findTimeIX(obj(1).time,delayedges);
delayt = delayt(1):delayt(2);
edges = [-0.4 -0.02]; % motion energy sort times
t = findTimeIX(obj(1).time,edges);

for isess = 2%1:numel(meta) % choice==2, ramp==1
    trialdat_zscored = permute(obj(isess).trialdat, [1 3 2]);
    trix{1} = (params(isess).trialid{cond2use(1)});
    trix{2} = (params(isess).trialid{cond2use(2)});


    dat = rez(isess).recon.null;
    W = cd_null(isess).cd_mode_orth(:,cdix);
    proj.null = tensorprod(dat,W,3,1);


    dat = rez(isess).recon.potent;
    W = cd_potent(isess).cd_mode_orth(:,cdix);
    proj.potent = tensorprod(dat,W,3,1);

    for j = 1:numel(trix)
        tempme{j} = me(isess).data(:,trix{j});
        nullts{j} = proj.null(:,trix{j});
        potentts{j} = proj.potent(:,trix{j});

        % sort trials by avg late delay motion energy
        [~,sortix] = sort(mean(tempme{j}(t(1):t(2),:),1),'ascend');

        tempme{j} = tempme{j}(delayt,sortix);
        nullts{j} = nullts{j}(delayt,sortix);
        potentts{j} = potentts{j}(delayt,sortix);

    end

    plt.me = cat(2,tempme{1},tempme{2});
    plt.null = cat(2,nullts{1},nullts{2});
    plt.potent = cat(2,potentts{1},potentts{2});

    plt.null = mySmooth(plt.null,nullsm,'reflect');
    plt.potent = mySmooth(plt.potent,potentsm,'reflect');

    % plt.null = zscore(plt.null);
    % plt.potent = zscore(plt.potent);

    
    f = figure;
    f.Position = [412   417   929   336];
    f.Renderer = 'painters';

    nTrials = size(plt.me,2);
    nCond = size(tempme{1},2);
    time_ = obj(1).time(delayt);

    ax = subplot(1,3,1); hold on;
    imagesc(time_,1:size(plt.me,2),plt.me'); colorbar; clim([5 80])
    ax = prettifyPlot(ax);
    colormap(ax,parula);
    plotEventTimes(ax,evtimes);
    line([time_(1) time_(end)], [nCond nCond],'Color','w','LineStyle','--');
    xlim([-2.39 0])
    ylim([1 numel(trix{1})]) % for cd ramp

    ax = subplot(1,3,2); hold on;
    imagesc(time_,1:size(plt.me,2),plt.potent'); %colorbar; %clim([-max(max(plt.potent))/1.5 max(max(plt.potent))/1.5])
    ax = prettifyPlot(ax);
    colormap(ax,cm);
    plotEventTimes(ax,evtimes);
    line([time_(1) time_(end)], [nCond nCond],'Color','k','LineStyle','-');
    xlim([-2.39 0])
    % ylim([1 numel(trix{1})]) % for cd ramp

    ax = subplot(1,3,3); hold on;
    imagesc(time_,1:size(plt.me,2),plt.null'); %colorbar; %clim([-10 50])%max(max(plt.null))])
    ax = prettifyPlot(ax);
    colormap(ax,cm);
    plotEventTimes(ax,evtimes);
    line([time_(1) time_(end)], [nCond nCond],'Color','k','LineStyle','--');
    xlim([-2.39 0])
    % ylim([1 numel(trix{1})]) % for cd ramp

    sgtitle([meta(isess).anm ' ' meta(isess).date])

    % break
    % pause
    
    
end



%%