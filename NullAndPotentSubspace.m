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
% 1) dynamically set number of dimensions per subspace based on some heuristic
%       could be upper bound dimensionality for each of the covariance matrices
%       (Cov-mov and Cov-nomov   or   Cov-resp and Cov-delay)
% 2) inpar.delay and inpar.responseOnly should be more robust. I just hardcoded
%       the time values right now
% 3) create inpar.standardize.method , which can be 'baseline' or 'entire
%       trial'
% 4) deal with dual probes for visualization and analysis after finding
%       subspaces
% 5) make subspace identification code very modular, so that in.method and
%       inpar.standardize and finding VE and projs and all that shares code across
%       each method of SID
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

meta = meta(1:2);

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'st'; % 'st' or 'ta' or '2pca' or 'regress'

% inpar.trials = 'all'; % specify 'all' or condition numbers
inpar.trials = 'all'; %2:5;

inpar.delayOnly = false; % if true, only use delay epoch for subspace estimation
inpar.responseOnly = false; % if true, only use response epoch for subspace estimation

% dimensionality (will soon change to dynamically set this)
inpar.nNullDim = 4;
inpar.nPotentDim = 4;

inpar.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)

inpar.estimateDimensionality = false; % if true, find upper-bound of dimensionality using parallel analysis
inpar.dimpth = 'results\Dimensionality'; % where to save parallel analysis results

%% LOAD DATA and PERFORM SUBSPACE ANALYSIS

all_alignment = struct();
all_ve.null = nan(size(meta));
all_ve.potent = nan(size(meta));
all_contrib = struct();
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth all_alignment all_ve all_contrib
    
    % load data
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta)) ' : ' meta(isess).anm ' ' meta(isess).date])
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
    
    % identify subspaces
    [in,out] = Subspace_SingleTrialElsayed(inpar,sessobj,me.move,sesspar,thismeta,params);
    
    % projections
    proj = ProjectDataToSubspace(in.data.zscored,out.Q);
    
    % reconstructions and variance explained
    for isub = 1:numel(inpar.subspace_names)
        thissub = inpar.subspace_names{isub};
        recon.(thissub) = ReconstructDataFromSubspace(proj.(thissub),out.Q.(thissub));
        [ve.unit.(thissub), ve.total.(thissub)] = SubspaceReconVarianceExplained(in.data.zscored, recon.(thissub));
    end

    % subspace alignment 
    mask = ve.unit.null > 0.05 & ve.unit.potent > 0.05; % null and potent subspace must explain 20% of data variance
    out.alignment.all = (ve.unit.null(mask) - ve.unit.potent(mask)) ./ (ve.unit.null(mask) + ve.unit.potent(mask));
    
    % tagged unit subspace alignment
    if numel(sessobj.psth)>1
        taggedids = tag.id.obj{1};
        taggedids = [taggedids  ; ( tag.id.obj{2}+taggedids(end) )];
    else
        taggedids = tag.id.obj{1};
    end
    out.alignment.pt = (ve.unit.null(taggedids) - ve.unit.potent(taggedids)) ./ (ve.unit.null(taggedids) + ve.unit.potent(taggedids));

    % subspace contribution
    contrib.null.all = sqrt(sum(out.Q.null.^2,2)) / size(out.Q.null,2) .* out.fr.null;
    totalcontrib = sum(contrib.null.all);
    contrib.null.all = contrib.null.all / totalcontrib;
    contrib.potent.all = sqrt(sum(out.Q.potent.^2,2)) / size(out.Q.potent,2) .* out.fr.potent;
    totalcontrib = sum(contrib.potent.all);
    contrib.potent.all = contrib.potent.all / totalcontrib;
    contrib.null.pt = contrib.null.all(taggedids);
    contrib.potent.pt = contrib.potent.all(taggedids);
    
    all_alignment(isess).all = out.alignment.all;
    all_alignment(isess).pt = out.alignment.pt;
    all_ve.null(isess) = ve.total.null;
    all_ve.potent(isess) = ve.total.potent;
    all_contrib(isess).null.all = contrib.null.all;
    all_contrib(isess).potent.all = contrib.potent.all;
    all_contrib(isess).null.pt = contrib.null.pt;
    all_contrib(isess).potent.pt = contrib.potent.pt;
    
end

%% plot alignment distributions
close all


cols = getColors;

f = figure;
f.Renderer = 'painters';
f.Position = [680   581   410   297];
ax = gca;
ax = prettifyAxis(ax);
hold on;

% plot all units alignment distribution
toplot = cell2mat({all_alignment.all}');

h = histogram(toplot,30,'edgecolor','none','Normalization','probability','Visible','off');
bars = h.Values;
binedges = h.BinEdges;

% x = find(binedges<0);
% b = bar(binedges(x),bars(x));
% b.BarWidth = 1;
% b.EdgeColor = 'none';
% b.FaceColor = cols.potent;
% 
% x = find(binedges>0);
% x = x(1:end-1);
% b = bar(binedges(x),bars(x));
% b.BarWidth = 1;
% b.EdgeColor = 'none';
% b.FaceColor = cols.null;

% plot pt units alignments
toplot = cell2mat({all_alignment.pt}');
h = histogram(toplot*-1.0,10,'edgecolor','none','Normalization','probability');
h.FaceColor = [0.3 0.3 0.3];
xl = plot([0 0],ax.YLim,'k--');
xl.LineWidth = 2;
p = patch([ax.XLim(1) 0 0 ax.XLim(1)], [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],'k');
p.FaceColor = cols.potent;
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';
p = patch([0 ax.XLim(2) ax.XLim(2) 0], [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)],'k');
p.FaceColor = cols.null;
p.FaceAlpha = 0.2;
p.EdgeColor = 'none';

% change stacking order
set(ax, 'Children', flipud(get(ax, 'Children')) )

xlabel('Subspace alignment')
ylabel('Probability')

%% plot subspace contribution
close all

toplot_null = [];
toplot_potent = [];
toplotpt_null = [];
toplotpt_potent = [];
for i = 1:numel(all_contrib)
    toplot_null = cat(1,toplot_null,all_contrib(i).null.all);
    toplot_potent = cat(1,toplot_potent,all_contrib(i).potent.all);
    toplotpt_null = cat(1,toplotpt_null,all_contrib(i).null.pt);
    toplotpt_potent = cat(1,toplotpt_potent,all_contrib(i).potent.pt);
end

toplot = {toplot_null,toplotpt_null,toplot_potent,toplotpt_potent}; % (null,nullpt, potent,potentpt)


cols = getColors;
c(1,:) = cols.null;
c(2,:) = cols.null./1.5;
c(3,:) = cols.potent;
c(4,:) = cols.potent./1.5;

f = figure;
f.Position = [943   363   292   335];
f.Renderer = "painters";
hold on;
ax = prettifyAxis(gca);
hold on;
xs = [1 2 4 5];
for i = 1:numel(xs)
    this = toplot{i};
    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(this)./i, 0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor',c(i,:), 'markeredgecolor','none')
end
ax.XTick = xs;
xticklabels({'null','nullpt','potent','potentpt'})
ylabel('Subspace contribution')

%% plot total variance explained by null and potent subspace
close all
clear xs_ ys_

cols = getColors;
cc(1,:) = cols.null;
cc(2,:) = cols.potent;
cc(3,:) = [ 0 0 0];

% plot nUnits
f = figure;
f.Renderer = 'painters';
f.Position = [726   414   271   349];
ax = prettifyAxis(gca); 
hold on;
xs = [1 2 3];
data = {all_ve.null,all_ve.potent};
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
    xx = simple_violin_scatter(xs(i)*ones(size(this)), tempthis, numel(meta), 0.2);
    % xx = xs(i)*ones(size(this));
    scatter(xx, this, 30,'filled', 'markerfacecolor',cc(i,:), 'markeredgecolor',[0.8 0.8 0.8])
    xs_(:,i) = xx;
    ys_(:,i) = this;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end

ax.XTick = xs;
xticklabels({'Null','Potent'})
ylabel('Frac. VE')
xlim([0.5,2.5])
ylim([0 0.8])



