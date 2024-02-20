function [params,obj] = processData(obj,params,prbnum, varargin)

% varargin only accepts one input right now -> params.behav_only. If 1,
% don't process clusters
if nargin > 3
    behav_only = varargin{1};
else
    behav_only = 0;
end

%% STANDARD ROUTINES
% find trials to use (only need to do this once)
if prbnum==1 || ~isfield(params,'trialid')
    params.trialid = findTrials(obj, params.condition);
    disp(' ')
    disp('--Trials Found')
    disp(' ')
end

% find clusters to use
params.cluid{prbnum} = findClusters({obj.clu{prbnum}(:).quality}', params.quality);
disp(' ')
disp('--Clusters Found')
disp(' ')

if behav_only
    return
end

% warp spikes post go cue according to median lick duration for each lick
if isfield(params,'timeWarp')
    if params.timeWarp
        [obj,obj.lick] = timeWarp(obj,params,prbnum);
        disp(' ')
        disp('--Time Warp Finished')
        disp(' ')
    end
end

% align spikes in every cluster to an event
obj = alignSpikes(obj,params,prbnum);
disp(' ')
disp('--Spikes Aligned')
disp(' ')

% get trial avg psth and single trial data
if ~isfield(params,'bctype') % boundary condition for smoothing
    params.bctype = 'none';
end
obj = getSeq(obj,params,prbnum);
disp(' ')
disp('--PSTHs and Single Trial Data Done')
disp(' ')


% reject units based on quality metrics
[obj,params.cluid{prbnum}] = QMRejection(obj,params.cluid{prbnum},params.qm,prbnum);
% [obj, params.cluid{prbnum}] = removeLowFRClusters(obj,params.cluid{prbnum},params.lowFR,prbnum);
disp(' ')
disp('--Rejected units based on quality metrics')
disp(' ')
% 
% % get mean fr and std dev for each cluster/trial type during presample (baseline fr)
% [obj.presampleFR{prbnum}, obj.presampleSigma{prbnum}] = baselineFR(obj,params,prbnum);
% disp(' ')
% disp('--Presample Statisitics Calculated')
% disp(' ')


end










