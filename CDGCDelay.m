clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% TODO
% - calculate firing rate change from dleay to go cue per unit
% - calculate CD gocue-delay

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

params.behav_only = 0;

params.useEarly = 0;

params.dt = 1/100;

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

%% LOAD DATA
% delay = [-2.145 -1.85];
delay = [-0.4 0];
resp = [0 0.4];

cond = [2 3];

allquality = [];
allregion = [];
zresp.mu = [];
zresp.sigma = [];
allcd = [];
proj.ptn = [];
proj.un = [];
for isess = 1:numel(meta)
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    

    activity = [];
    psth = [];
    for iprobe = 1:numel(sessobj.trialdat)
        activity = cat(2,activity,sessobj.trialdat{iprobe});
        allquality = cat(1,allquality,sesspar.quality{iprobe});
        allregion = cat(1,allregion, cellstr(repmat(thismeta.region{iprobe},numel(sesspar.quality{iprobe}),1)));
        psth = cat(2,psth,sessobj.psth{iprobe}(:,:,cond));
    end

    % calulate change in firing rate from delay to go cue
    idelay = findTimeIX(sessobj.time,delay,1);
    iresp = findTimeIX(sessobj.time,resp,1);


    delayActivity = squeeze(nanmean(activity(idelay,:,:),1)); % (neurons,trials)
    respActivity = squeeze(nanmean(activity(iresp,:,:),1));

    mu_delay = mean(delayActivity, 2); 
    std_delay = std(delayActivity, [], 2);
    z  = (respActivity - mu_delay) ./ std_delay;

    zresp.mu = cat(1,zresp.mu, nanmean(z,2));
    zresp.sigma = cat(1,zresp.sigma, nanstd(z,[],2));

    
    mu_resp = mean(respActivity, 2);
    std_resp = std(respActivity, [], 2);

    mu = cat(2,mu_delay,mu_resp);
    sd = cat(2,std_delay,std_resp);

    cd = ((mu(:,2)-mu(:,1)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
    cd = cd./sum(abs(cd)); % (ncells,1)

    allcd = cat(1,allcd,cd);

    ptn_psth = activity(:,tag.tagix_allprobes,:);
    ptn_proj = nanmean(tensorprod(ptn_psth,cd(tag.tagix_allprobes),2,1), 2);

    iun = ~ismember(1:size(psth,2),tag.tagix_allprobes);
    un_psth = activity(:,iun,:);
    un_proj = nanmean(tensorprod(un_psth,cd(iun),2,1), 2);

    proj.ptn = cat(3,proj.ptn,ptn_proj);
    proj.un = cat(3,proj.un,un_proj);


end

%%

close all

iptn = ismember(allquality,'tagged');
iun = ~ismember(allquality,'tagged');

toplot = {zresp.mu(iptn)  zresp.mu(iun)};

cols = getColors;
c(1,:) = cols.ptn;
c(2,:) = cols.un;

f = figure;
f.Position = [943   363   292   335];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;
xs = [1 2];
for i = 1:numel(xs)
    this = toplot{i};
    if i==2
        this = this - 0.2;
    end
    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, min(numel(this),1000), 0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor',c(i,:), 'markeredgecolor','none')
end
ax.XTick = xs;
xticklabels({'PTN','Other'})
ylabel('z-scores from delay to response')


%%

f = figure;
f.Position = [943   363   292   335];
f.Renderer = "painters";
ax = prettifyAxis(gca);
hold on;

cols = getColors;

sel = squeeze(proj.ptn(:,2,:) - proj.ptn(:,1,:));
[m,h] = mean_CI(sel, 0.95);
shadedErrorBar(sessobj.time,m,h,{'Color',cols.ptn,'LineStyle','-','LineWidth',2},0.3,ax);

sel = squeeze(proj.un(:,2,:) - proj.un(:,1,:));
[m,h] = mean_CI(sel, 0.95);
shadedErrorBar(sessobj.time,m,h,{'Color',cols.un,'LineStyle','-','LineWidth',2},0.3,ax);



















