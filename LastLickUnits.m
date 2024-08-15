clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc


%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
% params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};
params.qm.quality = {'tagged'};

params.alignEvent = 'goCue';
params.tmin = -1;
params.tmax = 7;

params.condition = {};
params.condition(1)         = {'R&hit'};                   % right hit trials
params.condition(end+1)     = {'L&hit'};                   % left hit trials
params.condition(end+1)     = {'R&miss'};                  % right hit trials
params.condition(end+1)     = {'L&miss'};                  % left hit trials

params.condLabel = {};
params.condLabel{end+1} = 'rhit';
params.condLabel{end+1} = 'lhit';
params.condLabel{end+1} = 'rmiss';
params.condLabel{end+1} = 'lmiss';

params.behav_only = 0;

params.dt = 1/20;

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

%% subset meta (TODO)

% meta = subsetMetaByParams(meta,params);


%% LOAD DATA

for isess = 1:numel(meta)
    disp(' ')
    disp(['Session ' num2str(isess) '/' num2str(numel(meta))])
    disp(' ')
    [sessobj,sesspar] = loadSessionData(meta(isess),params);

    tag(isess) = getTagFromObj(sessobj,sesspar,meta(isess));
    
end

%% SAVE TAG

% % % fpth = 'C:\Users\munib\Documents\Economo-Lab\data\tagged';
% % % % fn = 'AllTagged_GoCue_20240817.mat';
% % % fn = 'AllTagged_LastLick_20240817.mat';
% % % % fn = 'AllTagged_FirstLick_20240817.mat';
% % % save(fullfile(fpth,fn),'tag','params','-v7.3')


%% LAST LICK UNITS
% somethings very wrong here. not getting right trials, lick times....

close all

cond = 1:2;
sm = 21;

for isess = 1:numel(tag)

    clear last_lick_time gocue_times all_lick_times lick_times

    % get all lick times
    % subtract go cue time from lick times
    % remove any negative lick times
    % find the last lick time
    last_lick_time = nan(tag(isess).bp.Ntrials,1);
    for itrial = 1:tag(isess).bp.Ntrials
        lickR = tag(isess).bp.ev.lickR{itrial};
        lickL = tag(isess).bp.ev.lickL{itrial};

        all_lick_times = sort([lickR lickL]);
        if isempty(all_lick_times)
            last_lick_time(itrial) = nan;
            continue
        end

        gocue_times(itrial) = tag(isess).bp.ev.goCue(itrial);
        
        all_lick_times = all_lick_times - gocue_times(itrial);
        all_lick_times(all_lick_times<0) = [];

        
        if isempty(all_lick_times)
            last_lick_time(itrial) = nan;
            continue
        else
            last_lick_time(itrial) = max(all_lick_times);
        end
    end

    % for each tag(isess)ged unit, for right and left hit trials
    for iprobe = 1:numel(tag(isess).nTag)
        for iunit = 1:tag(isess).nTag(iprobe)
            f = figure;
            f.Position = [680   141   389   737];
            f.Renderer = 'painters';
            t = tiledlayout('flow');
            thisunit = squeeze(tag(isess).trialdat{iprobe}(:,iunit,:)); % (time,trials
            for icond = 1:numel(cond)
                ax = nexttile;
                hold on;
                thiscond = cond(icond);
                thistrials = sesspar.trialid{thiscond};
                % sort trials by last_lick_time
                this_lick_times = last_lick_time(thistrials);
                [sorted_lick_times,sorted_trials] = sort(this_lick_times);
                toplot = thisunit(:,sorted_trials);
                toplot = mySmooth(toplot,sm,'reflect');
                imagesc(tag(isess).time,1:numel(sorted_trials),toplot')
                plot([0 0],ax.YLim,'w--','LineWidth',2)
                scatter(sorted_lick_times,1:numel(sorted_trials),10,'MarkerFaceColor','w','MarkerEdgeColor','none')
                % for i = 1:numel(sorted_
                title(params.condLabel{thiscond})
                ylim([1,numel(sorted_trials)])
                xlim([params.tmin,params.tmax])
                colormap(jet)
            end
        end
    end
end





















