clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'manopt')));
addpath(genpath(fullfile(utilspth,'subspace')));

clc

%% TODO
% 1) dimensionality of each subspace...
% 2) inpar.delay and inpar.responseOnly should be more robust. I just hardcoded
%       the time values right now
% 3) deal with dual probes for visualization and analysis after finding
%       subspaces

% something's wrong with projected variance explained for ta and two-pca

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single'};

params.behav_only = 0;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
resultspth = fullfile(utilspth,'results','subspaceID');
meta = [];

% meta = allSessionMeta(meta,datapth);
meta1 = ALM_SessionMeta(meta,datapth);
meta2 = tjM1_SessionMeta(meta,datapth);
meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth);  % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'pls'; % 'st' or 'ta' or 'two-pca' or 'regress' or 'pls'
inpar.regress.regularize = 'lasso'; % 'lasso' or 'ridge' or 'none (TODO)'
inpar.regress.pca = 0;  % use pcs instead of neurons for regression
inpar.regress.cv = 1; % cross validate to find best lasso or ridge param

% inpar.trials = 'all'; % specify 'all' or condition numbers
inpar.trials = 'all'; % 2:5;

inpar.delayOnly = false; % if true, only use delay epoch for subspace estimation
inpar.responseOnly = false; % if true, only use response epoch for subspace estimation

% dimensionality (will soon change to dynamically set this)
inpar.nNullDim = 8;
inpar.nPotentDim = 8;

inpar.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)

inpar.estimateDimensionality = false; % if true, find upper-bound of dimensionality using parallel analysis
inpar.dimpth = 'results\Dimensionality'; % where to save parallel analysis results

inpar.standardize = 'baseline'; % 'baseline' or 'trial' - which time points to compute zscore stats over

inpar.saveResults = 0;

nIters = 50;

%% LOAD DATA and PERFORM SUBSPACE ANALYSIS

clear W
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth all_alignment all_ve all_contrib resultspth nIters W

    % load data
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta)) ' : ' meta(isess).anm ' ' meta(isess).date])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
    allquality = sesspar.quality{:};
    qual_mask = ismember(allquality,'single');

    for iter = 1:nIters
        if mod(iter,10)==0
            disp(['Iteration: ' num2str(iter) '/' num2str(nIters)])
        end
        % perform subspace identification
        [in,out] = SubspaceEngine(inpar,sessobj,sesspar,me);

        % % collect weights of tagged units and all units
        w = sqrt(sum(out.W.^2,2));

        W.full.all{isess}(:,iter) = w(qual_mask);
        W.full.pt{isess}(:,iter) = w(tag.tagix_allprobes);

        wp = sqrt(sum(out.Q.potent.^2,2));
        wn = sqrt(sum(out.Q.null.^2,2));

        W.potent.all{isess}(:,iter) = wp(qual_mask);
        W.null.all{isess}(:,iter) = wn(qual_mask);

        W.potent.pt{isess}(:,iter) = wp(tag.tagix_allprobes);
        W.null.pt{isess}(:,iter) = wn(tag.tagix_allprobes);
    end
end

%% plot subspace contribution (potent and null)
close all
clear allW ptW

cols = getColors;

allW.potent = vertcat(W.potent.all{:});
ptW.potent = vertcat(W.potent.pt{:});
allW.null = vertcat(W.null.all{:});
ptW.null = vertcat(W.null.pt{:});

f = figure;
f.Position = [943   414   732   284];
f.Renderer = "painters";
ax = prettifyAxis(subplot(1,2,1));
hold on;

histogram(allW.potent,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.un)
histogram(ptW.potent,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.ptn)
title('Potent')
xlabel('norm(Wpotent)');
ylabel('Probability Density');
leg = legend('Single units', 'PTNs'); leg.Box = 'off';

ax = prettifyAxis(subplot(1,2,2));
hold on;

histogram(allW.null,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.un)
histogram(ptW.null,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.ptn)
title('Null')
xlabel('norm(Wnull)');

%% plot subspace contribution (just the W matrix)
close all
clear allW ptW

cols = getColors;

allW = vertcat(W.full.all{:});
% allW = normalize(allW(:),'range',[0 1]);
ptW = vertcat(W.full.pt{:});
% ptW = normalize(ptW(:),'range',[0 1]);

allW = allW(:);
allWix = 1:numel(allW);
ptW = ptW(:);
ptWix = 1:numel(ptW);

WW = cat(1,allW,ptW);
WW = normalize(WW,'range',[0 1]);
allW = WW(allWix);
ptW = WW(ptWix+allWix(end));



f = figure;
f.Position = [943   405   385   293];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;

histogram(allW,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.un)
histogram(ptW,8,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.ptn)
xlabel('Neuron contribution');
ylabel('Probability Density');
leg = legend('Single units', 'PTNs'); leg.Box = 'off';
title('Lasso regression','fontsize',11.5)

% % Create inset for range [0.6, 1]
% inset_ax = axes('Position', [0.5, 0.38, 0.33, 0.3]);  % Adjust [x, y, width, height] as needed
% inset_ax = prettifyAxis(inset_ax);
% hold(inset_ax, 'on');
% 
% % Filter data for the specified range
% 
% % Plot histograms on the inset axis
% hall = histogram(inset_ax, allW, 'Normalization', 'pdf', 'EdgeAlpha', 0, 'FaceColor', cols.un);
% hpt = histogram(inset_ax, ptW, 'Normalization', 'pdf', 'EdgeAlpha', 0, 'FaceColor', cols.ptn);
% 
% % Set inset axis limits and labels
% xlim([0.55 1.05])
% ylim([0 0.5])


%% plot relative subspace contribution 
close all
clear allW ptW

cols = getColors;

allW.potent = vertcat(W.potent.all{:});
ptW.potent = vertcat(W.potent.pt{:});
allW.null = vertcat(W.null.all{:});
ptW.null = vertcat(W.null.pt{:});

allW.rel = (allW.null - allW.potent) ./ (allW.null + allW.potent);
ptW.rel = (ptW.null - ptW.potent) ./ (ptW.null + ptW.potent);


f = figure;
f.Position = [943   414   732   284];
f.Renderer = "painters";
ax = prettifyAxis(subplot(1,2,1));
hold on;

histogram(allW.rel,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.un)
histogram(ptW.rel,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.ptn)
title('Potent')
xlabel('norm(Wpotent)');
ylabel('Probability Density');
leg = legend('Single units', 'PTNs'); leg.Box = 'off';

ax = prettifyAxis(subplot(1,2,2));
hold on;

histogram(allW.null,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.un)
histogram(ptW.null,'Normalization','pdf','EdgeAlpha',0,'FaceColor',cols.ptn)
title('Null')
xlabel('norm(Wnull)');
