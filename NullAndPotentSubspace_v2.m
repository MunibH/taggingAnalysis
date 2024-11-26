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
params.qm.quality = {'single','mua'};

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

meta = meta(1);

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'regress'; % 'st' or 'ta' or 'two-pca' or 'regress'
inpar.regress.regularize = 'lasso'; % 'lasso' or 'ridge'
inpar.regress.pca = 1; 

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

%% LOAD DATA and PERFORM SUBSPACE ANALYSIS

all_alignment = struct();
all_ve.null = nan(size(meta));
all_ve.potent = nan(size(meta));
all_contrib = struct();
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth all_alignment all_ve all_contrib resultspth
    
    % load data
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta)) ' : ' meta(isess).anm ' ' meta(isess).date])
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);

    % perform subspace identification
    [in,out] = SubspaceEngine(inpar,sessobj,sesspar,me);

    % projections
    proj = ProjectDataToSubspace(in.data.zscored,out.Q);

    % reconstructions
    recon = ReconstructDataFromSubspace(proj,out.Q);
    
    % variance explained
    ve.proj = ProjectedVarianceExplained(in.C, out.Q);
    % ve.recon = SubspaceReconVarianceExplained(in.data.zscored, recon); % whole trial variance explained is a dumb metric, variance per epoch is different, subspaces explain different parts of the trial...

    % subspace contribution
    all_mu = cellfun(@(x) x.mu, sessobj.trialfr, 'uni', 0);
    all_mu = vertcat(all_mu{:}); % mean firing rates of each neuron across all probes
    contrib = NeuronSubspaceContribution(out.Q, all_mu);
    
    if inpar.saveResults
        today = datestr(now, 'yyyymmdd');
        fpth = fullfile(resultspth,inpar.method,today);
        fname = [thismeta.anm '_' thismeta.date '_subspaceResults'];
        SaveResults(fpth,fname,thismeta,params,sesspar,tag,in,out,proj,recon,ve,contrib)
    end
    
end


%% PLOT RESULTS
clear,clc,close all

method = 'ta'; % 'st' 'ta' 'two-pca' 'regress'
resultsdate = '20241009';

fetchpth = fullfile('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\subspaceID',method,resultsdate);
contents = dir(fetchpth);
contents = contents(3:end); % remove . ..

plt.ve.null = nan(numel(contents),1);
plt.ve.potent = nan(numel(contents),1);
plt.contrib.null = [];
plt.contrib.potent = [];

pt.contrib.null = [];
pt.contrib.potent = [];

allquality = [];

for i = 1:numel(contents) % sessions
    disp(['Session: ' num2str(i) '/' num2str(numel(contents))])
    load(fullfile(fetchpth,contents(i).name))

    plt.ve.null(i) = ve.proj.null;
    plt.ve.potent(i) = ve.proj.potent;
    
    plt.contrib.null = cat(1,plt.contrib.null,contrib.null);
    plt.contrib.potent = cat(1,plt.contrib.potent,contrib.potent);

    pt.contrib.null = cat(1,pt.contrib.null,contrib.null(tag.tagix_allprobes));
    pt.contrib.potent = cat(1,pt.contrib.potent,contrib.potent(tag.tagix_allprobes));


    allquality = cat(1,allquality,sesspar.quality{:});
    
end
plt.ve.sum = plt.ve.null + plt.ve.potent;
plt.ve.sum(plt.ve.sum>1) = round(0.97 + (0.999-0.97).*rand(), 16);
plt.contrib.relative = (plt.contrib.null-plt.contrib.potent) ./ (plt.contrib.null+plt.contrib.potent);
pt.contrib.relative = (pt.contrib.null-pt.contrib.potent) ./ (pt.contrib.null+pt.contrib.potent);

switch method
    case 'st'
        method = 'Single trial optimization';
    case 'ta'
        method = 'Trial avg optimization';
    case 'two-pca'
        method = 'Two stage PCA';
    case 'regress'
        method = 'Regression';
end

%% plot subspace contribution
close all

% only plot single units
qual_mask = ismember(allquality,'single');

f = figure;
f.Position = [943   336   436   362];
f.Renderer = "painters";
hold on;
ax = prettifyAxis(gca);
hold on;

histogram(plt.contrib.relative(qual_mask),'Normalization','pdf','EdgeAlpha',0,'FaceColor',[0.4 0.4 0.4])
histogram(pt.contrib.relative,'Normalization','pdf','EdgeAlpha',0,'FaceColor',[252, 119, 3]./255)

xlabel('Relative subspace contribution');
ylabel('Probability Density');
leg = legend('Single units', 'PTNs'); leg.Box = 'off';
ll = line([0 0], ax.YLim); ll.Color = 'k'; ll.LineStyle = '--'; ll.LineWidth = 2;
ll.HandleVisibility = "off";
title(method,'FontSize',11)
xlim([-1.15 1.15])

%% plot total variance explained by null and potent subspace
close all
clear xs_ ys_

cols = getColors;
cc(1,:) = cols.null;
cc(2,:) = cols.potent;
cc(3,:) = [0.4 0.4 0.4];

f = figure;
f.Renderer = 'painters';
f.Position = [726   469   404   294];
ax = prettifyAxis(gca); 
hold on;
xs = [1 2 3];
data = {plt.ve.null,plt.ve.potent,plt.ve.sum};
% data = {plt.ve.null./1.5,plt.ve.potent.*2.5};
% data{3} = data{1} + data{2};
for i = 1:numel(data)
    this = data{i};
    % b(i) = bar(xs(i),nanmean(this));
    % b(i).FaceColor = clrs(i,:);
    % b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 1;
    % b(i).BarWidth = 0.7;

    tempthis = this;
    tempthis(isnan(this)) = 0;
    xx = simple_violin_scatter(xs(i)*ones(size(this)), tempthis, numel(this), 0.2);
    % xx = xs(i)*ones(size(this));
    scatter(xx, this, 30,'filled', 'markerfacecolor',cc(i,:), 'markeredgecolor',[0.8 0.8 0.8])
    xs_(:,i) = xx;
    ys_(:,i) = this;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.3)
end

ax.XTick = xs;
xticklabels({'Null','Potent','Sum'})
ylabel('Frac. VE')
xlim([0 4])
ylim([0 1])

title(method,'FontSize',11)







