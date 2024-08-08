function [params,obj] = processData(obj,meta,params,prbnum,edges)


%% STANDARD ROUTINES

% find clusters to use
[maskQuality,allquality] = findClustersByQuality(obj.clu{prbnum},obj.metrics{prbnum}, params.qm.quality);
clumask = maskQuality;
disp('~~~~~~~~~~~ Clusters Found ~~~~~~~~~~~')

% reject units based on quality metrics
if params.qm.perform
    maskQM = QMRejection(obj,params.qm,prbnum);
    clumask = clumask&maskQM;
    % [obj, params.cluid{prbnum}] = removeLowFRClusters(obj,params.cluid,params.lowFR,prbnum);
    disp('~~~~~~~~~~~ Rejected units based on quality metrics ~~~~~~~~~~~')
end

if params.removeTagOverlap
    maskTagOverlap = rejectTagOverlap(obj,params.tagOverlapThresh, prbnum);
    clumask = clumask&maskTagOverlap;
    disp('~~~~~~~~~~~ Rejected units based on overlap with tagged units ~~~~~~~~~~~')
end


params.clumask = clumask;
params.cluid = find(params.clumask);

% retrieve which shank, channel, and quality for each unit in params.cluid
[params.channel,params.shank] = getShankAndChannels(obj,prbnum,params.cluid);
params.quality = allquality(params.clumask);
params.region = repmat({meta.region{prbnum}},numel(params.cluid),1);

if params.behav_only
    return
end

% warp spikes post go cue according to median lick duration for each lick
if params.timeWarp
    obj = timeWarp2(obj,params,prbnum);
    disp('~~~~~~~~~~~ Time Warp Finished ~~~~~~~~~~~')
end

% align spikes in every cluster to an event
obj = alignSpikes(obj,params,prbnum);
disp(['~~~~~~~~~~~ Spikes Aligned to ' params.alignEvent ' ~~~~~~~~~~~'])

% get trial avg psth and single trial data
if ~isfield(params,'bctype') % boundary condition for smoothing
    params.bctype = 'none';
end
obj = getSeq(obj,params,prbnum,edges);
disp('~~~~~~~~~~~ PSTHs and Single Trial Data Done ~~~~~~~~~~~')




% 
% % get mean fr and std dev for each cluster/trial type during presample (baseline fr)
[obj.baseline.mu, obj.baseline.sigma] = baselineFR(obj,params,prbnum);
disp('~~~~~~~~~~~ Presample Statisitics Calculated ~~~~~~~~~~~')


end










