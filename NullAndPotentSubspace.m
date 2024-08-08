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

meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'st'; % 'st' or 'ta' or '2pca' or 'regress'

inpar.trials = 'all'; % specify 'all' or condition numbers
% inpar.trials = 2:5;

inpar.delayOnly = false; % if true, only use delay epoch for subspace estimation
inpar.responseOnly = false; % if true, only use response epoch for subspace estimation

% dimensionality (will soon change to dynamically set this)
inpar.nNullDim = 4;
inpar.nPotentDim = 4;

inpar.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)


%% LOAD DATA and PERFORM SUBSPACE ANALYSIS

all_alignment = struct();
all_ve.null = nan(size(meta));
all_ve.potent = nan(size(meta));
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth all_alignment all_ve
    
    % load data
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta))])
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
    
    all_alignment(isess).all = out.alignment.all;
    all_alignment(isess).pt = out.alignment.pt;
    all_ve.null(isess) = ve.total.null;
    all_ve.potent(isess) = ve.total.potent;
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
h = histogram(toplot*-1,10,'edgecolor','none','Normalization','probability');
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

    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(meta), 0.2);
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
% ylim([0 0.8])



