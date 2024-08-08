function [obj,params] = loadSessionData(meta,params,varargin)

if nargin > 2
    onlyObj = varargin{1};
else
    onlyObj = false;
end

%% load data obj
disp(' ')
disp(['~~~~~~~~~~~ LOADING DATA OBJ: ' meta.anm ' ' meta.date ' ~~~~~~~~~~~'])
disp(' ')
load(meta.datapth); % loads 'obj'
if onlyObj
    params = nan;
    return
end

%% delete tagging trials from obj if tagging session
if params.remove_tag
    disp("~~~~~~~~~~~ Omitting tagging trials ~~~~~~~~~~~")
    obj = deleteTaggingTrials(obj);
end

%% trials by condition
disp('~~~~~~~~~~~ Finding trials ~~~~~~~~~~~')
params.trialid = findTrials(obj, params.condition);

%% time
edges = params.tmin:params.dt:params.tmax;
obj.time = edges + params.dt/2;
obj.time = obj.time(1:end-1)';

%% add fake probe for testing
% meta.probe = [1 2];
% obj.clu{2} = obj.clu{1}(1:20);


%% process data per specified probe
params.probe = meta.probe;
params.region = meta.region;

for prbix = 1:numel(params.probe)
    disp(['~~~~~~~~~~~ Processing data for Probe ' num2str(params.probe(prbix)) ', ' meta.region{prbix} ' ~~~~~~~~~~~'])

    prbnum = params.probe(prbix);
    
    % find clusters
    % time warp
    % align spikes
    % firing rates
    % unit rejection via quality metrics

    [probeparams{prbix},probeobj{prbix}] = processData(obj,meta,params,prbnum,edges);
end

%% 
if numel(params.probe) == 1

    params = probeparams{1};
    obj = probeobj{1};

    params.clumask = {params.clumask};
    params.cluid = {params.cluid};
    params.shank = {params.shank};
    params.channel = {params.channel};
    params.quality = {params.quality};
    params.region = {params.region};
    
    if params.behav_only
    else
        obj.psth = {obj.psth};
        obj.trialdat = {obj.trialdat};
        obj.baseline = {obj.baseline};
    end


elseif numel(params.probe) == 2

    params.clumask = {probeparams{1}.clumask, probeparams{2}.clumask};
    params.cluid = {probeparams{1}.cluid, probeparams{2}.cluid};
    params.shank = {probeparams{1}.shank, probeparams{2}.shank};
    params.channel = {probeparams{1}.channel, probeparams{2}.channel};
    params.quality = {probeparams{1}.quality, probeparams{2}.quality};
    params.region = {probeparams{1}.region, probeparams{2}.region};

    if params.behav_only
    else
        obj.psth = {probeobj{1}.psth, probeobj{2}.psth};
        obj.trialdat = {probeobj{1}.trialdat, probeobj{2}.trialdat};
        obj.bp = probeobj{1}.bp;
        obj.baseline = {probeobj{1}.baseline, probeobj{2}.baseline};
    end

else 
    error('more than dual probe processing not yet implemented')
end

%% get event times
params.eventTimes = getEventTimes(obj.bp.ev,params.events,params.alignEvent);

disp(' ')
disp('~~~~~~~~~~~ DATA LOADED AND PROCESSED ~~~~~~~~~~~')
disp(' ')



