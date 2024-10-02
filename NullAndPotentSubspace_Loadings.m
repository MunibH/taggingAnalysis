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
% meta1 = ALM_SessionMeta(meta,datapth);
% meta2 = tjM1_SessionMeta(meta,datapth);
% meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth);  % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

meta = meta(1);

%% PARAMETERS

inpar.subspace_names = {'null','potent'};

inpar.method = 'st'; % 'st' or 'ta' or '2pca' or 'regress'

% inpar.trials = 'all'; % specify 'all' or condition numbers
inpar.trials = 2:5;

inpar.delayOnly = false; % if true, only use delay epoch for subspace estimation
inpar.responseOnly = false; % if true, only use response epoch for subspace estimation

% dimensionality (will soon change to dynamically set this)
inpar.nNullDim = 4;
inpar.nPotentDim = 4;

inpar.alpha = 1e100; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)

inpar.estimateDimensionality = false; % if true, find upper-bound of dimensionality using parallel analysis
inpar.saveDimensionality = false;
inpar.dimpth = 'results\Dimensionality'; % where to save parallel analysis results

nBoots = 1;

%% LOAD DATA and PERFORM SUBSPACE ANALYSIS
clear contrib Q thismeta sessobj sesspar tag me

sess_boot_ct = 1;
for isess = 1:numel(meta)
    clearvars -except isess meta inpar params datapth utilspth iboot nBoots sess_boot_ct Q contrib
    % load data
    disp(['Session: ' num2str(isess) '/' num2str(numel(meta))])
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
    for iboot = 1:nBoots
        disp(['Iteration: ' num2str(iboot) '/' num2str(nBoots)])
        clearvars -except isess meta inpar params datapth utilspth iboot nBoots sess_boot_ct sessobj sesspar tag me thismeta Q contrib

        bootobj = sessobj;
        bootpar = sesspar;
        inpar.nProbes = numel(bootobj.trialdat);

        % identify subspaces
        [in,out] = Subspace_SingleTrialElsayed(inpar,bootobj,me.move,bootpar,thismeta);
        
        % save loadings
        Q.null(:,:,sess_boot_ct) = out.Q.null;
        Q.potent(:,:,sess_boot_ct) = out.Q.potent;

        % subspace contribution
        % con = sqrt(sum(out.Q.null.^2,2)) / size(out.Q.null,2) .* out.fr.null;
        con = sqrt(sum(out.Q.null.^2,2));% .* out.fr.null;
        totalcontrib = sum(con);
        contrib.null(:,sess_boot_ct) = con / totalcontrib;
        % con = sqrt(sum(out.Q.potent.^2,2)) / size(out.Q.potent,2) .* out.fr.potent;
        con = sqrt(sum(out.Q.potent.^2,2));% .* out.fr.potent;
        totalcontrib = sum(con);
        contrib.potent(:,sess_boot_ct) = con / totalcontrib;

        sess_boot_ct = sess_boot_ct + 1;
    end
end

%% are the loadings gaussian? (yes, centered around 0)
close all

cols = getColors;

n = Q.null(:);
p = Q.potent(:);

f = figure;
% f.Position = [18         458        1894         420];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
histogram(n,100,'FaceColor',cols.null,'EdgeColor','none','Normalization','probability','FaceAlpha',0.6)
histogram(p,100,'FaceColor',cols.potent,'EdgeColor','none','Normalization','probability','FaceAlpha',0.4)

xlabel('Loading')
ylabel('Probability')

%% what about subspace contribution? (half-normal or maybe exponential, mean close to 0)

close all

cols = getColors;

n = contrib.null(:);
p = contrib.potent(:);

f = figure;
% f.Position = [18         458        1894         420];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
histogram(n,50,'FaceColor',cols.null,'EdgeColor','none','Normalization','probability','FaceAlpha',0.6)
histogram(p,50,'FaceColor',cols.potent,'EdgeColor','none','Normalization','probability','FaceAlpha',0.4)

xlabel('Contribution')
ylabel('Probability')

%% how does dimensionality effect these results? (it doesn't)


%% plot variability in contribution (l2 norm of loadings without FR)

close all

cols = getColors;
nUnits = sum(cell2mat(cellfun(@numel,sesspar.cluid,'uni',0)));

f = figure;
f.Position = [18         458        1894         420];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;

temp = contrib.null;
boxplot(temp','PlotStyle','compact','Colors',cols.null,'Whisker',5,...
    'Positions',1:nUnits)
% set(gca,'XTickLabel',{' '})

temp = contrib.potent;
boxplot(temp','PlotStyle','compact','Colors',cols.potent,'Whisker',5,...
    'Positions',(1:nUnits)+0.2)
% set(gca,'XTickLabel',{' '})
xlabel('Unit')
ylabel('Contribution (w/out FR)')


%% plot some intersting units
close all

unitnums = [1 8 13 17 28 47 55 61 69 90 91 93 96 115];

q = cat(1,sesspar.quality{:});

cond = [2 3];

psth = cat(2,sessobj.psth{:});

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

f = figure;
% f.Position = [18         458        1894         420];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
for i = 1:numel(unitnums)% size(psth,2)
    cla(ax)
    u = unitnums(i);
    % u = i;
    for icond = 1:numel(cond)
        temp = psth(:,u,cond(icond));
        plot(sessobj.time,temp,'Color',c(icond,:),'linewidth',2)
    end
    xlim([-2.1,params.tmax])
    plotEventTimes(ax,sesspar.eventTimes,'k',false);
    title(['Null: ' num2str(round(mean(contrib.null(u,:)),3)) ...
        ', Potent: ' num2str(round(mean(contrib.potent(u,:)),3)) ... 
        ', ' q{u}], 'fontsize',10.5)
    pause

end

















