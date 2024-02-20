%% NP alignment and PSTHs

% find L/R selective cells per session
edges = [-0.8 0];

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');
cols = getColors;

close all
psth.null = [];
psth.potent = [];
psth.no = [];
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

    cond2use = [5 6]; %[8 9];
    psth.null = cat(2,psth.null,obj(sessix).psth(:,clus(nmask),cond2use));
    psth.potent = cat(2,psth.potent,obj(sessix).psth(:,clus(pmask),cond2use));
    psth.no = cat(2,psth.no,obj(sessix).psth(:,clus(nomask),cond2use));

    
    

    % % PLOT PSTHS AND ALIGNMENT
    % f = figure;
    % f.Position = [425         173        1061         750];
    % t = tiledlayout('flow');
    % for i = 1:numel(alignment)
    %     ax = nexttile;
    %     hold on;
    %     ax = prettifyPlot(ax);
    % 
    %     thisclu = clus(i);
    % 
    %     plot(obj(sessix).time,obj(sessix).psth(:,thisclu,5),'Color',cols.rhit,'linewidth',2)
    %     plot(obj(sessix).time,obj(sessix).psth(:,thisclu,6),'Color',cols.lhit,'linewidth',2)
    %     title([num2str(round(alignment(i),2)) ' ' meta(sessix).anm ' ' meta(sessix).date])
    %     % title(['Alignment = ' num2str(round(alignment(i),2))])
    %     xlabel('Time from go cue (s)')
    %     ylabel('Firing rate')
    %     plotEventTimes(ax,evtimes)
    %     xlim([-2.3,2])
    % 
    % end

end


%% pca for each category (null,potent,no aligned)

fns = fieldnames(psth);
for i = 1:numel(fns)
    this = fns{i};

    dat = psth.(this);
    dims = size(dat);

    temp = reshape(permute(dat,[1 3 2]),dims(1)*dims(3),dims(2));

    temp = zscore(temp);

    pc = pca(temp,'NumComponents',3);

    p = temp*pc;
    % p = zscore(temp) * pc;

    proj.(this) = reshape(p,dims(1),dims(3),3);

    % proj.(this) = tensorprod(dat,pc,2,1); %(time,cond,pc)

end

% plot



evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');

cols = getColors;
clrs(1,:) = cols.null;
clrs(2,:) = cols.potent;
clrs(3,:) = [0.4 0.4 0.4];

    
f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
ax = subplot(3,3,1);
ax = prettifyPlot(ax);
hold on;
this = 'null';
mode = 1;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);
% ylabel('PC1')

ax = subplot(3,3,2);
ax = prettifyPlot(ax);
hold on;
this = 'potent';
mode = 1;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);

ax = subplot(3,3,3);
ax = prettifyPlot(ax);
hold on;
this = 'no';
mode = 1;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);

%
ax = subplot(3,3,4);
ax = prettifyPlot(ax);
hold on;
this = 'null';
mode = 2;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);
% ylabel('PC2')

ax = subplot(3,3,5);
ax = prettifyPlot(ax);
hold on;
this = 'potent';
mode = 2;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);

ax = subplot(3,3,6);
ax = prettifyPlot(ax);
hold on;
this = 'no';
mode = 2;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);

%
ax = subplot(3,3,7);
ax = prettifyPlot(ax);
hold on;
this = 'null';
mode = 3;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);
% ylabel('PC3')

ax = subplot(3,3,8);
ax = prettifyPlot(ax);
hold on;
this = 'potent';
mode = 3;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);

ax = subplot(3,3,9);
ax = prettifyPlot(ax);
hold on;
this = 'no';
mode = 3;
toplot = proj.(this)(:,:,mode);
plot(obj(1).time,toplot(:,1),'Color',cols.rhit,'linestyle','-','linewidth',2)
plot(obj(1).time,toplot(:,2),'Color',cols.lhit,'linestyle','-','linewidth',2)
plotEventTimes(ax,evtimes)
xlim([-2.39 2]);


