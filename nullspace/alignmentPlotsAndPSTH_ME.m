%% NP alignment and PSTHs

% find L/R selective cells per session
edges = [-0.8 0];

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');
cols = getColors;

close all

mee = [-0.8,0.5];
meix = findTimeIX(obj(1).time,mee);
meix = meix(1):meix(2);
percentiles = [30 70];

for sessix = 1:numel(obj)
    cond2use = [5 6]; %[8 9];
    cluix = findSelectiveCells(obj(sessix),params(sessix),edges,cond2use);

    cond2use = {'hit|miss'};
    clear trialdat W proj r2
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix);

    trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:);

    W.null = rez(sessix).Qnull;
    W.potent = rez(sessix).Qpotent;
    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.recon.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)

        for k = 1:numel(clus) % for each cell
            thisclu = clus(k);

            % calculcate variance explained by CD choice
            orig = trialdat.full(:,:,thisclu); % (time,trials)

            fr = mean(mean(orig)); % subspace contribution method
            % weight = norm(W.(fns{j})(k));
            % r2.(fns{j}){sessix}(k) = fr*weight;

            recon = trialdat.recon.(fns{j})(:,:,thisclu); % (time,trials) % ve by recon method
            mdl = fitlm(orig(:),recon(:));
            r2.(fns{j})(k) = mdl.Rsquared.Ordinary;
        end
    end
    alignment = (r2.null - r2.potent) ./ (r2.null + r2.potent);

    cond2use = [5 6]; %[8 9];
    conds = params(1).condition(cond2use);
    trix = findTrials(obj(sessix), conds);
    
    data.m = cellfun(@(x) me(sessix).data(:,x), trix, 'uni',0);
    data.n = cellfun(@(x) obj(sessix).trialdat(:,:,x), trix, 'uni',0);

    % find high/low ME trials
    medat = cellfun(@(x) mean(x(meix,:),1), data.m, 'uni', 0);
    p = cellfun(@(x) prctile(x,percentiles), medat, 'uni', 0);

    lowtrix{1} = medat{1} <= p{1}(1); % right low me trials
    lowtrix{2} = medat{2} <= p{2}(1); % left low me trials

    hightrix{1} = medat{1} >= p{1}(2); % right high me trials
    hightrix{2} = medat{2} >= p{2}(2); % left high me trials


    % PLOT PSTHS AND ALIGNMENT
    f = figure;
    f.Position = [425         173        1061         750];
    t = tiledlayout('flow');
    for i = 1:numel(alignment)
        if (alignment(i) >= -0.7) && (alignment(i) <=0.7)
            continue
        end
        ax = nexttile;
        hold on;
        ax = prettifyPlot(ax);

        thisclu = clus(i);
        

        this.low.right.mu = mean(squeeze(data.n{1}(:,thisclu,lowtrix{1})),2);
        this.low.left.mu = mean(squeeze(data.n{2}(:,thisclu,lowtrix{2})),2);
        this.high.right.mu = mean(squeeze(data.n{1}(:,thisclu,hightrix{1})),2);
        this.high.left.mu = mean(squeeze(data.n{2}(:,thisclu,hightrix{2})),2);

        % this.low.right.sem = std(squeeze(data.n{1}(:,thisclu,lowtrix{1})),[],2)./sqrt(numel(lowtrix{1}));
        % this.low.left.sem = std(squeeze(data.n{2}(:,thisclu,lowtrix{2})),[],2)./sqrt(numel(lowtrix{2}));
        % this.high.right.sem = std(squeeze(data.n{1}(:,thisclu,hightrix{1})),[],2)./sqrt(numel(hightrix{1}));
        % this.high.left.sem = std(squeeze(data.n{2}(:,thisclu,hightrix{2})),[],2)./sqrt(numel(hightrix{2}));
        
        this.low.right.sem = getCI(squeeze(data.n{1}(:,thisclu,lowtrix{1})),1);
        this.low.left.sem = getCI(squeeze(data.n{2}(:,thisclu,lowtrix{2})),1);
        this.high.right.sem = getCI(squeeze(data.n{1}(:,thisclu,hightrix{1})),1);
        this.high.left.sem = getCI(squeeze(data.n{2}(:,thisclu,hightrix{2})),1);

        shadedErrorBar(obj(sessix).time,this.low.right.mu,this.low.right.sem,...
            {'Color',[103, 103, 235]./255,'LineWidth',2},0.3,ax);
        shadedErrorBar(obj(sessix).time,this.low.left.mu,this.low.left.sem,...
            {'Color',[240, 129, 144]./255,'LineWidth',2},0.3,ax);
        shadedErrorBar(obj(sessix).time,this.high.right.mu,this.high.right.sem,...
            {'Color',[11, 11, 99]./255,'LineWidth',2},0.3,ax);
        shadedErrorBar(obj(sessix).time,this.high.left.mu,this.high.left.sem,...
            {'Color',[99, 15, 26]./255,'LineWidth',2},0.3,ax);

        % plot(obj(sessix).time,this.low.right,'Color',[103, 103, 235]./255,'linewidth',2)
        % plot(obj(sessix).time,this.low.left,'Color',[240, 129, 144]./255,'linewidth',2)
        % plot(obj(sessix).time,this.high.right,'Color',[11, 11, 99]./255,'linewidth',2)
        % plot(obj(sessix).time,this.high.left,'Color',[99, 15, 26]./255,'linewidth',2)
        title([num2str(round(alignment(i),2)) ' ' meta(sessix).anm ' ' meta(sessix).date])
        % title(['Alignment = ' num2str(round(alignment(i),2))])
        xlabel('Time from go cue (s)')
        ylabel('Firing rate')
        plotEventTimes(ax,evtimes)
        xlim([-2.3,2])

    end



end





