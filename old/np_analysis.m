clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - tprime opto out
% - time warp
% - handle tagged units post loading data

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

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% subset meta (TODO)

% meta = subsetMetaByParams(meta,params);


%% LOAD DATA

thismeta = meta(1);

[sessobj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(sessobj,sesspar,thismeta);

me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
kin = getKinematics(sessobj, me, sesspar);

%% NULL SPACE

nNPDims = 8;% nNPDims as last arg if you want to change nDims per subspace [4,6,10,13]

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
disp('finding null and potent spaces')
% -- input data
trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
% trialdat_zscored = permute(obj.trialdat, [1 3 2]);

% -- null and potent spaces
cond2use = [2,3,4,5]; % right hit, left hit, right miss, left miss, 2afc
cond2proj = 1:5; % ~early versions
nullalltime = 0; % use all time points to estimate null space if 1
onlyAW = 0; % only use AW trials
delayOnly = 0; % only use delay period
responseOnly = 0;
rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), ...
    params(sessix), cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly,nNPDims); % nNPDims as last arg if you want to change nDims per subspace


disp('DONE')

%% projections

dat.null = tensorprod(trialdat_zscored,rez.Qnull,3,1);
dat.potent = tensorprod(trialdat_zscored,rez.Qpotent,3,1);

%% reconstructions

recon.null = tensorprod(dat.null,rez.Qnull,3,2);
recon.potent = tensorprod(dat.potent,rez.Qpotent,3,2);

%% plot reconstructions
close all

f = figure;
f.Position = [680   646   767   232];
f.Renderer = 'painters';
ax1 = prettifyAxis(subplot(1,2,1));
hold(ax1,'on');
ax2 = prettifyAxis(subplot(1,2,2));
hold(ax2,'on');

cols = getColors;
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

plt = 1;
if plt

    cond2plot = [2,3];
    for i = rez.N.dims(3)-tag.nTag+1:rez.N.dims(3)
        cla(ax1)
        cla(ax2)
        c(1,:) = cols.rhit;
        c(2,:) = cols.lhit;
        for j = 1:numel(cond2plot)
            trix = params.trialid{cond2plot(j)};
            this = squeeze(mean(trialdat_zscored(:,trix,i),2));
            plot(ax1,obj.time,this,'color',c(j,:),'LineWidth',2)
            plot(ax2,obj.time,this,'color',c(j,:),'LineWidth',2)
        end
        c(1,:) = cols.rhit_aw*1.2;
        c(2,:) = cols.lhit_aw*1.2;
        c(c>1) = 1;
        for j = 1:numel(cond2plot)
            trix = params.trialid{cond2plot(j)};
            thisn = squeeze(mean(recon.null(:,trix,i),2));
            thisp = squeeze(mean(recon.potent(:,trix,i),2));
            plot(ax1,obj.time,thisn,'color',c(j,:),'LineWidth',3)
            plot(ax2,obj.time,thisp,'color',c(j,:),'LineWidth',3)
        end


        plotEventTimes(ax1,params.eventTimes,'k',false)
        plotEventTimes(ax2,params.eventTimes,'k',false)
        xlim(ax1,[params.tmin,params.tmax])
        xlim(ax2,[params.tmin,params.tmax])
        title(ax1,'Null')
        title(ax2,'Potent')

        pause

    end
end

%% ve

fns = {'null','potent'};
for j = 1:numel(fns)
    % single trials neural activity reconstructed from n/p
    rc = recon.(fns{j});

    for k = 1:size(rc,3) % for each cell
        orig = trialdat_zscored(:,:,k); % (time,trials)

        pred = rc(:,:,k); % (time,trials) % ve by recon method
        mdl = fitlm(orig(:),pred(:));
        r2.(fns{j})(k) = mdl.Rsquared.Ordinary;
    end
end

%%

alignment = (r2.null - r2.potent) ./ (r2.null + r2.potent);

f = figure;
f.Position = [680   581   410   297];
ax = gca;
ax = prettifyAxis(ax);
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

tagalignment = alignment(end-tag.nTag+1:end);
mx = ax.YLim(2) - 4;

scatter(tagalignment,ones(size(tagalignment))*mx, 50,'k','filled','markeredgecolor','w')
xlabel('subspace alignment')
ylabel('unit count')


%% save alignment data

save(['tagalignment_' meta.anm '_' meta.date '.mat'],'tagalignment')



