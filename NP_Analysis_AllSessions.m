clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'nullspace')));

clc

% % TODO
% - selectivity
% - correlations and xcorrelations with movements
% - kinematic decoding like in our paper
% - np analysis
% - - need to automate finding ME threshold

%% Parameters

datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
allmeta = [];
allmeta = allSessionMeta(allmeta,datapth);

% specify changes to params here
chgparams.alignEvent = 'firstLick';

% how many NP dims to find
nNPDims = 8;

% save results from null and potent analysis
saveResults = 0;
savePath = fullfile(utilspth,'results','NullPotent');
savePath = fullfile(savePath,datestr(now, 'yyyymmdd_HHMMSS'));

% save figures
saveReconFig = 1;
saveAlignmentFig = 0;
saveCorrelationFig = 0;

%% analyze single sessions, save results

for isess = 1:numel(allmeta)

    clearvars -except allmeta datapth chgparams nNPDims ...
        saveResults savePath saveReconFig saveAlignFig saveCorrelationFig isess

    % get default params, set any changes
    params = defaultParams();
    if exist('defparams','var')
        chgparamsf = fieldnames(chgparams);
        for i = 1:numel(chgparamsf)
            this = chgparamsf{i};
            params.(this) = chgparams.(this);
        end
    end

    % get current session
    meta = allmeta(isess);

    params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

    % LOAD DATA

    [obj,params] = loadSessionData(meta,params);
    me = loadMotionEnergy(obj, meta, params, datapth);
    kin = getKinematics(obj, me, params);

    % TAGGED UNIT META
    tag.nTag = numel(obj.tag);
    tag.cluid = [obj.tag(:).cluid]; % where tagged units are in obj.clu
    tag.cluid = find(ismember(params.cluid,tag.cluid))'; % where in params.cluid, trialdat, psth

    % NULL SPACE
    % -----------------------------------------------------------------------
    % -- Curate Input Data --
    % zscore single trial neural data (time*trials,neurons), for all trials
    % -- Calculate null and potent spaces --
    % null space from quiet time points
    % potent space from moving time points
    % -----------------------------------------------------------------------
    disp('finding null and potent spaces')
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj);
    % trialdat_zscored = permute(obj.trialdat, [1 3 2]);

    % -- null and potent spaces
    cond2use = [2,3,4,5]; % right hit, left hit, right miss, left miss, 2afc
    cond2proj = 1:5; % ~early versions
    nullalltime = 0; % use all time points to estimate null space if 1
    onlyAW = 0; % only use AW trials
    delayOnly = 0; % only use delay period
    responseOnly = 0;
    rez = singleTrial_elsayed_np(trialdat_zscored, obj, me, ...
        params, cond2use, cond2proj, nullalltime, onlyAW, delayOnly, responseOnly,nNPDims); % nNPDims as last arg if you want to change nDims per subspace


    disp('DONE')

    % projections
    rez.trialdat_zscored = trialdat_zscored;
    rez.proj.null = tensorprod(trialdat_zscored,rez.Qnull,3,1);
    rez.proj.potent = tensorprod(trialdat_zscored,rez.Qpotent,3,1);

    % reconstructions
    rez.recon.null = tensorprod(rez.proj.null,rez.Qnull,3,2);
    rez.recon.potent = tensorprod(rez.proj.potent,rez.Qpotent,3,2);

    % VE

    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        rc = rez.recon.(fns{j});

        for k = 1:size(rc,3) % for each cell
            orig = trialdat_zscored(:,:,k); % (time,trials)

            pred = rc(:,:,k); % (time,trials) % ve by recon method
            temp = corrcoef(orig(:),pred(:));
            r2.(fns{j})(k) = temp(1,2).^2;
        end
    end

    % ALIGNMENT
    alignment = (r2.null - r2.potent) ./ (r2.null + r2.potent);
    tagalignment = alignment(end-tag.nTag+1:end);

    % save alignment data
    rez.alignment.nontag = alignment;
    rez.alignment.tag = tagalignment;
    rez.params = params;
    rez.meta = meta;
    rez.tag = tag;

    if saveResults
        saveFn = [meta.anm '_' meta.date '_NP.mat'];
        SaveResults(savePath,saveFn,'rez',rez);
    end


end


%% load all sessions results
clearvars -except savePath

% override savePath if you want to load specific results
savePath = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\NullPotent\20240417_145401';

% get all data of interest
% just going to plot alignment data
filelist = dir(fullfile(savePath,'*.mat'));  %get list of files and folders in any subfolder containing token
filelist = filelist(~[filelist.isdir]);  %remove folders from list

% data I want from each results file
align.nontag = [];
align.tag = [];
regions = cell(1);
for isess = 1:numel(filelist)
    clear a
    a = load(fullfile(filelist(isess).folder,filelist(isess).name)); % loads rez
    rez(isess) = a.rez;

    align.nontag = [align.nontag a.rez.alignment.nontag];
    align.tag = [align.tag a.rez.alignment.tag];

    ntag = a.rez.tag.nTag;
    for i = 1:ntag
        regions = cat(1,regions,a.rez.meta.region);
    end

end
clear a
regions = regions(2:end);

%% plot alignment

close all

cols = getColors;

f = figure;
f.Position = [680   581   410   297];
ax = gca;
ax = prettifyAxis(ax,4);
hold on;

h = histogram(align.nontag,30,'edgecolor','none','Normalization','count','Visible','off');
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

mx = ax.YLim(2) - 10;

scatter(align.tag,ones(size(align.tag))*mx, 200, '|', 'MarkerEdgeColor','k')
ln = line([0 0], ax.YLim);
ln.LineStyle = '--';
ln.Color = 'k';
ln.LineWidth = 2;
xlabel('Subspace alignment')
ylabel('Unit count')

if saveAlignmentFig
    mysavefig(f,savePath,'alignment')
end


%% plot reconstructions

% plot reconstructions
close all

cols = getColors;

cond2plot = [2 , 3];

for isess = 1:numel(rez)
    clear thisrez
    thisrez = rez(isess);

    time = getTimeVector(thisrez.params.tmin,thisrez.params.tmax,thisrez.params.dt);

    c(1,:) = cols.rhit;
    c(2,:) = cols.lhit;

    taggedUnits = thisrez.N.dims(3)-thisrez.tag.nTag+1:thisrez.N.dims(3);

    for itag = 1:numel(taggedUnits) % for every tagged unit
        i = taggedUnits(itag);

        f = figure;
        f.Position = [680   646   767   232];
        f.Renderer = 'painters';
        ax1 = prettifyAxis(subplot(1,2,1));
        hold(ax1,'on');
        ax2 = prettifyAxis(subplot(1,2,2));
        hold(ax2,'on');

        cla(ax1)
        cla(ax2)
        c(1,:) = cols.rhit;
        c(2,:) = cols.lhit;
        for j = 1:numel(cond2plot)
            trix = thisrez.params.trialid{cond2plot(j)};
            this = squeeze(mean(thisrez.trialdat_zscored(:,trix,i),2));
            plot(ax1,time,this,'color',c(j,:),'LineWidth',2)
            plot(ax2,time,this,'color',c(j,:),'LineWidth',2)
        end
        c(1,:) = cols.rhit_aw*1.2;
        c(2,:) = cols.lhit_aw*1.2;
        c(c>1) = 1;
        for j = 1:numel(cond2plot)
            trix = thisrez.params.trialid{cond2plot(j)};
            thisn = squeeze(mean(thisrez.recon.null(:,trix,i),2));
            thisp = squeeze(mean(thisrez.recon.potent(:,trix,i),2));
            plot(ax1,time,thisn,'color',c(j,:),'LineWidth',3)
            plot(ax2,time,thisp,'color',c(j,:),'LineWidth',3)
        end


        plotEventTimes(ax1,thisrez.params.eventTimes,'k',true)
        plotEventTimes(ax2,thisrez.params.eventTimes,'k',true)
        xlim(ax1,[thisrez.params.tmin,thisrez.params.tmax])
        xlim(ax2,[thisrez.params.tmin,thisrez.params.tmax])

        algn = round(thisrez.alignment.tag(itag),2);
        if algn < 0
            titleax = ax2;
        else
            titleax = ax1;
        end
        title(titleax,['Alignment: ' num2str(algn)],'fontsize',12)

        ylabel(ax1,'Null (zscore)')
        ylabel(ax2,'Potent (zscore)')
        xlabel(ax1,'Time from go cue (s)')

        if saveReconFig
            unitid = thisrez.params.cluid(thisrez.tag.cluid(itag));
            thisfn = ['NPReconPSTH_' thisrez.meta.anm '_' thisrez.meta.date '_Unit' num2str(unitid) ];
            mysavefig(f,savePath,thisfn,1)
        end

    end

end



%% how many null/potent aligned units in each region

alm_mask = contains(regions,{'ALM'});
tjm1_mc_mask = contains(regions,{'tjM1','MC'});

N.alm.potent = sum(alm_mask' & (align.tag<0));
N.tjm1.potent = sum(tjm1_mc_mask' & (align.tag<0));
N.alm.null = sum(alm_mask' & (align.tag>0));
N.tjm1.null = sum(tjm1_mc_mask' & (align.tag>0));


%%
% I eventually want to answer the question:
% are the strongly potent aligned units overrpresented

% there's a confound though
% the alignment distribution is strongly potent skewed
% and there may be many non-tagged units in there that are correlated with
% the activity of tagged units. Those non-tagged units may be 'missed' PT
% cells... thus making the distribution representative of PT cells and the
% comparison would be nonsignificant

% for now, just gonna calculate the correlation between each unit and each
% tagged unit

cond2use = [2,3];
psths = nan(rez(1).N.dims(1),19,numel(cond2use)); % (time,tagunits,conds)

cc = nan(2000,19); % (max num of nontag units, num tag units)
ii = 1;
for isess = 1:numel(rez)
    clear thisrez
    thisrez = rez(isess);

    taggedUnits = thisrez.N.dims(3)-thisrez.tag.nTag+1:thisrez.N.dims(3);
    nNonTagged = thisrez.N.dims(3) - numel(taggedUnits);

    for itag = 1:numel(taggedUnits) % for every tagged unit
        i = taggedUnits(itag);

        tagdat = thisrez.trialdat_zscored(:,:,i); % (time,trials)

        for j = 1:nNonTagged
            nontagdat = thisrez.trialdat_zscored(:,:,j); % (time,trials)
            tempcc = corrcoef(tagdat(:),nontagdat(:));
            cc(j,ii) = tempcc(1,2);
        end


        % get psths
        for j = 1:numel(cond2use)
            trix = thisrez.params.trialid{cond2use(j)};
            psths(:,ii,j) = squeeze(mean(thisrez.trialdat_zscored(:,trix,i),2));
        end

        ii = ii + 1;
    end
end


%% % plot correlations
close all

f = figure;
f.Renderer = 'painters';
f.Position = [5         490        1911         315];
ax = prettifyAxis(gca,0.5);
hold on;

cols = linspecer(size(cc,2));

xs = 1:size(cc,2);
for i = 1:numel(xs)
    this = cc(:,i);
    % mu = nanmean(this);
    % b(i) = bar(xs(i),mu);
    % b(i).FaceColor = cols(i,:);
    % b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 1;
    % b(i).BarWidth = 0.7;

    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(rez), 0.8);
    scatter(xx, this, 20, 'markerfacecolor',cols(i,:), 'markeredgecolor','none')

    thispsth = mySmooth(squeeze(psths(:,i,:)),15,'none',0);
    thispsth = normalize(thispsth,'range',[-0.9,-0.45]);
    xs_ = linspace(xs(i)-0.5,xs(i)+0.4,size(thispsth,1));
    plot(xs_,thispsth(:,1),'b','LineWidth',1)
    plot(xs_,thispsth(:,2),'r','LineWidth',1)


end

ax.XTick = xs;
ylabel('Correlation on single trials')
xlabel('Tagged unit number')
ylim([ax.YLim(1),1])
ylim([-1,1])



if saveCorrelationFig
    thisfn = 'CorrTagNonTag_AllTaggedUnits';
    mysavefig(f,savePath,thisfn,0)
end


%% plot trial averaged data correlations
% 
% 
% % now correlate psths, for every tagged unit with all other units across
% % sessions and within session
% 
% 
% tagpsths = psths; % from previous section (time,units,conds)
% nontagpsths = nan(size(psths,1),5000,size(psths,3)); % (time,num nontag units,conds)
% 
% cond2use = [2,3];
% ii = 1;
% for isess = 1:numel(rez)
%     clear thisrez
%     thisrez = rez(isess);
% 
%     taggedUnits = thisrez.N.dims(3)-thisrez.tag.nTag+1:thisrez.N.dims(3);
%     nNonTagged = thisrez.N.dims(3) - numel(taggedUnits);
% 
% 
%     % get psths
%     for i = 1:nNonTagged
%         for j = 1:numel(cond2use)
%             trix = thisrez.params.trialid{cond2use(j)};
%             nontagpsths(:,ii,j) = squeeze(mean(thisrez.trialdat_zscored(:,trix,i),2));
%         end
%         ii = ii + 1;
%     end
% end
% 
% nontagpsths = nontagpsths(:,1:ii,:); % remove empty
% 
% % correlate psths
% nTag = size(tagpsths,2);
% nNonTag = size(nontagpsths,2);
% cc = nan(nNonTag,nTag);
% for itag = 1:nTag
%     for inontag = 1:nNonTag
%         tag_ = squeeze(tagpsths(:,itag,:));
%         nontag_ = squeeze(nontagpsths(:,inontag,:));
% 
%         tempcc = corrcoef(tag_(:),nontag_(:));
%         cc(inontag,itag) = tempcc(1,2);
%     end
% end
% 
% close all
% 
% f = figure;
% f.Renderer = 'painters';
% f.Position = [5         490        1911         315];
% ax = prettifyAxis(gca,0.5);
% hold on;
% 
% cols = linspecer(size(cc,2));
% 
% xs = 1:size(cc,2);
% for i = 1:numel(xs)
%     this = cc(:,i);
% 
%     xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(rez), 0.8);
%     scatter(xx, this, 5, 'markerfacecolor',cols(i,:), 'markeredgecolor','none')
% 
%     thispsth = mySmooth(squeeze(psths(:,i,:)),15,'none',0);
%     thispsth = normalize(thispsth,'range',[-1.5,-0.9]);
%     xs_ = linspace(xs(i)-0.5,xs(i)+0.4,size(thispsth,1));
%     plot(xs_,thispsth(:,1),'b','LineWidth',1)
%     plot(xs_,thispsth(:,2),'r','LineWidth',1)
% 
% end
% 
% ax.XTick = xs;
% ylabel('Correlation (PSTH)')
% xlabel('Tagged unit number')
% ylim([ax.YLim(1),1])
% ylim([-1.6,1])
% 
% 






