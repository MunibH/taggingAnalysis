%% find cells to use

clear cluix q
% find L/R selective cells per session
edges = [-0.8 0];
cond2use = [5 6]; %[8 9];
for i = 1:numel(obj)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
end

%% alignment NP
% reconstruct single cell activity from n/p

clear r2

cond2use = {'hit|miss'};
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});


    % trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    % trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);

    trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:);

    % rez(sessix).recon.null = tensorprod(rez(sessix).null_trialdat,rez(sessix).Qnull,3,2);
    % rez(sessix).recon.potent = tensorprod(rez(sessix).potent_trialdat,rez(sessix).Qpotent,3,2);

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
            r2.(fns{j}){sessix}(k) = mdl.Rsquared.Ordinary;
        end
    end
end


% plot

alln = [];
allp = [];
for sessix = 1:(numel(r2.null))
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
end

alignment = (alln - allp) ./ (allp + alln);

cols = getColors;

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;

x = find(binedges>0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;

% savepth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3\figNP\figs\NPalignment';
% fn = ['NPalignment_' num2str(nDims) ];
% saveas(f,fullfile(savepth,fn),'svg')

%% CD Choice alignment
% reconstruct single cell activity from n/p cd choice

clear r2

cdix = 1;
cond2use = {'hit|miss'};
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});


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
            r2.(fns{j}){sessix}(k) = mdl.Rsquared.Ordinary;
        end
    end
end


% plot

alln = [];
allp = [];
for sessix = 1:(numel(r2.null))
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
end

alignment = (alln - allp) ./ (allp + alln);

cols = getColors;

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;

x = find(binedges>0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;
% 
% savepth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3\figNP\figs\choiceAlignment';
% fn = ['ChoiceAlignment_' num2str(nDims) ];
% saveas(f,fullfile(savepth,fn),'svg')

%% CD Ramp alignment
% reconstruct single cell activity from n/p cd choice

clear r2

cdix = 3;
cond2use = {'hit|miss'};
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});


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
            r2.(fns{j}){sessix}(k) = mdl.Rsquared.Ordinary;
        end
    end
end


% plot

alln = [];
allp = [];
for sessix = 1:(numel(r2.null))
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
end

alignment = (alln - allp) ./ (allp + alln);

cols = getColors;

f = figure;
f.Position = [644   483   338   231];
ax = gca;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,30,'edgecolor','none','Normalization','count','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

x = find(binedges<0);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.potent;

x = find(binedges>0);
x = x(1:end-1);
b = bar(binedges(x),bars(x));
b.BarWidth = 1;
b.EdgeColor = 'none';
b.FaceColor = cols.null;


% savepth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3\figNP\figs\rampAlignment';
% fn = ['RampAlignment_' num2str(nDims) ];
% saveas(f,fullfile(savepth,fn),'svg')
