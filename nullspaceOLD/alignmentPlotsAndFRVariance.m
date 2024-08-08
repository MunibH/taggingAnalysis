%% NP alignment and PSTHs

% find L/R selective cells per session
edges = [-0.8 0];

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');
cols = getColors;

close all
psth.null = [];
psth.potent = [];
psth.no = [];

delaye = [-0.8 0];
delayix = findTimeIX(obj(1).time,delaye);
goe = [0 0.8];
goix = findTimeIX(obj(1).time,goe);

v.null.delay = [];
v.null.go = [];
v.potent.delay = [];
v.potent.go = [];
v.no.delay = [];
v.no.go = [];
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

    nmask = alignment > 0.6;
    pmask = alignment < -0.6;
    nomask = (alignment > -0.05) & (alignment < 0.05);

    data.null = obj(sessix).trialdat(:,clus(nmask),trix);
    data.potent = obj(sessix).trialdat(:,clus(pmask),trix);
    data.no = obj(sessix).trialdat(:,clus(nomask),trix);

    fns = fieldnames(data);
    for i = 1:numel(fns)
        this = fns{i};
        temp = data.(this);

        v.(this).delay = cat(2,v.(this).delay,mean(var(temp(delayix(1):delayix(2),:,:),[],1),3)); %(var over time,mean over trials, for each neuron)
        v.(this).go = cat(2,v.(this).go,mean(var(temp(goix(1):goix(2),:,:),[],1),3)); 

    end


end


%%


lab = {'null','potent','no'};

f=figure;
f.Position = [680   715   318   263];
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;

cols = getColors;
clrs(1,:) = cols.null;
clrs(2,:) = cols.potent;
clrs(3,:) = [0.4 0.4 0.4];

nCols = 3;
xs = 1:3;
fns = fieldnames(v);
for i = 1:3
    this = fns{i};
    temp = cat(1,v.(this).delay,v.(this).go)';
    bar(xs(i),nanmean(temp,1));
    % break
    % b(i).FaceColor = clrs(i,:);
    % b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 0.7;
    % b(i).BarWidth = 0.7;
    % 
    % vio = getXCoordsForViolin(this, []);
    % xx = (vio.jitter.*vio.jitterStrength./2);
    % 
    % xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(meta), 0.5);
    % 
    % % xs(i)*ones(numel(meta),1) + xx
    % scatter(xx, this, 20, 'markerfacecolor',clrs(i,:), 'markeredgecolor','k', 'linewidth',1)
    % % vs(i) = scatter(randn(numel(anms),1) * 0.1 + xs(i)*ones(size(perf(:,i))),perf(:,i),20,'MarkerFaceColor',cols{i},...
    % %     'MarkerEdgeColor','k','LineWidth',1);%,'XJitter','randn','XJitterWidth',0.25);
    % errorbar(b(i).XEndPoints,nanmean(this),getCI(this,1)./sqrt(numel(meta)),'LineStyle','none','Color','k','LineWidth',1)
end

ax.XTick = xs;
xticklabels(lab)
ax = gca;
ax.FontSize = 12;











