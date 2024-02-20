%% NP alignment and Depth

cols = getColors;

% find L/R selective cells per session
edges = [-0.8 0];

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');

close all
f = figure;
ax = gca;
ax = prettifyPlot(ax);
hold on;
dat.a.null = [];
dat.a.potent = [];
dat.d.null = [];
dat.d.potent = [];
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

    clear d
    % get depth
    depth = getDepth(meta(sessix),obj(sessix));
    % subset by clus of interest
    if ~iscell(depth)
        if isnan(depth)
            disp([meta(sessix).anm ' ' meta(sessix).date])
            continue
        end
    end
    if numel(meta(sessix).probe) > 1
        for iprobe = 1:numel(meta(sessix).probe)
            d{iprobe} = depth{iprobe}(params(sessix).cluid{iprobe});
        end
        d = cat(1,d{:});
    else
        d = depth(params(sessix).cluid); % subset based on params.cluid
    end

    d = d(clus); % subset based on selective cells
    
    pmask = alignment <= 0;
    nmask = alignment > 0;

    plot(alignment(pmask),d(pmask),'.','MarkerSize',10,'Color',cols.potent)
    plot(alignment(nmask),d(nmask),'.','MarkerSize',10,'Color',cols.null)

    dat.a.null = cat(1,dat.a.null,alignment(nmask)');
    dat.a.potent = cat(1,dat.a.potent,alignment(pmask)');
    dat.d.null = cat(1,dat.d.null,d(nmask));
    dat.d.potent = cat(1,dat.d.potent,d(pmask));

end
xline(0,'k--')
xlabel('NP Alignment')
ylabel('Depth (um)')

%% compare distributions

thresh = 0.6;
nmask = dat.a.null >=thresh;
pmask = dat.a.potent <= -thresh;

ndepth = dat.d.null(nmask);
pdepth = dat.d.potent(pmask);

[h,p] = ttest2(ndepth,pdepth)
% [p,h] = ranksum(ndepth,pdepth)


%% NP CDchoice alignment and Depth

cols = getColors;

% find L/R selective cells per session
edges = [-0.8 0];

cdix = 1;

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');

close all
f = figure;
ax = gca;
ax = prettifyPlot(ax);
hold on;
dat.a.null = [];
dat.a.potent = [];
dat.d.null = [];
dat.d.potent = [];
for sessix = 1:numel(obj)
    cond2use = [5 6]; %[8 9];
    cluix = findSelectiveCells(obj(sessix),params(sessix),edges,cond2use);

    cond2use = {'hit|miss'};
    clear trialdat W proj r2
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix);

    trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    % trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);
    

    W.null = cd_null(sessix).cd_mode_orth(:,cdix);
    W.potent = cd_potent(sessix).cd_mode_orth(:,cdix);
    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)
        % project onto Wchoice
        proj.(fns{j}) = tensorprod(trialdat.(fns{j}),W.(fns{j}),3,1);
        % reconstruct n/p reconstructions from CD choice proj
        trialdat.recon.(fns{j}) = tensorprod(proj.(fns{j}),W.(fns{j}),3,2);

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

    clear d
    % get depth
    depth = getDepth(meta(sessix),obj(sessix));
    % subset by clus of interest
    if ~iscell(depth)
        if isnan(depth)
            disp([meta(sessix).anm ' ' meta(sessix).date])
            continue
        end
    end
    if numel(meta(sessix).probe) > 1
        for iprobe = 1:numel(meta(sessix).probe)
            d{iprobe} = depth{iprobe}(params(sessix).cluid{iprobe});
        end
        d = cat(1,d{:});
    else
        d = depth(params(sessix).cluid); % subset based on params.cluid
    end

    d = d(clus); % subset based on selective cells
    
    pmask = alignment <= 0;
    nmask = alignment > 0;

    plot(alignment(pmask),d(pmask),'.','MarkerSize',10,'Color',cols.potent)
    plot(alignment(nmask),d(nmask),'.','MarkerSize',10,'Color',cols.null)

    dat.a.null = cat(1,dat.a.null,alignment(nmask)');
    dat.a.potent = cat(1,dat.a.potent,alignment(pmask)');
    dat.d.null = cat(1,dat.d.null,d(nmask));
    dat.d.potent = cat(1,dat.d.potent,d(pmask));

end
xline(0,'k--')
xlabel('CDchoice Alignment')
ylabel('Depth (um)')

%% compare distributions

thresh = 0.6;
nmask = dat.a.null >=thresh;
pmask = dat.a.potent <= -thresh;

ndepth = dat.d.null(nmask);
pdepth = dat.d.potent(pmask);

[h,p] = ttest2(ndepth,pdepth)
% [p,h] = ranksum(ndepth,pdepth)












