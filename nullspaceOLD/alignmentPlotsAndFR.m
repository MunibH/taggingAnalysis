%% NP alignment and Firing Rate

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
dat.f.null = [];
dat.f.potent = [];
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

    clear fr
    % get firing rate
    epoch = 'response'; % sample, delay , or response
    fr = getFiringRate(obj(sessix),epoch);
    fr = fr(clus); % subset based on selective cells
    
    pmask = alignment <= 0;
    nmask = alignment > 0;

    plot(alignment(pmask),fr(pmask),'.','MarkerSize',10,'Color',cols.potent)
    plot(alignment(nmask),fr(nmask),'.','MarkerSize',10,'Color',cols.null)

    dat.a.null = cat(1,dat.a.null,alignment(nmask)');
    dat.a.potent = cat(1,dat.a.potent,alignment(pmask)');
    dat.f.null = cat(1,dat.f.null,fr(nmask));
    dat.f.potent = cat(1,dat.f.potent,fr(pmask));

end
xline(0,'k--')
xlabel('NP Alignment')
ylabel('Firing Rate')

%% compare distributions

thresh = 0.6;
nmask = dat.a.null >=thresh;
pmask = dat.a.potent <= -thresh;

nfr = dat.f.null(nmask);
pfr = dat.f.potent(pmask);

[h,p] = ttest2(nfr,pfr)
% [p,h] = ranksum(ndepth,pdepth)

%% NP CDchoice alignment and Firing Rate

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
dat.f.null = [];
dat.f.potent = [];
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

    clear fr
    % get firing rate
    epoch = 'delay'; % sample, delay , or response
    fr = getFiringRate(obj(sessix),epoch);
    fr = fr(clus); % subset based on selective cells
    
    pmask = alignment <= 0;
    nmask = alignment > 0;

    plot(alignment(pmask),fr(pmask),'.','MarkerSize',10,'Color',cols.potent)
    plot(alignment(nmask),fr(nmask),'.','MarkerSize',10,'Color',cols.null)

    dat.a.null = cat(1,dat.a.null,alignment(nmask)');
    dat.a.potent = cat(1,dat.a.potent,alignment(pmask)');
    dat.f.null = cat(1,dat.f.null,fr(nmask));
    dat.f.potent = cat(1,dat.f.potent,fr(pmask));


end
xline(0,'k--')
xlabel('CDchoice Alignment')
ylabel('Firing Rate')

%% compare distributions

thresh = 0.6;
nmask = dat.a.null >=thresh;
pmask = dat.a.potent <= -thresh;

nfr = dat.f.null(nmask);
pfr = dat.f.potent(pmask);

[h,p] = ttest2(nfr,pfr)
% [p,h] = ranksum(ndepth,pdepth)



