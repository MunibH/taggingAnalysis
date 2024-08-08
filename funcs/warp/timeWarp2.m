function [obj] = timeWarp2(obj,params,prbnum)

% use jaw phase to warp each lick cycle post go cue
view = 1; % side cam
coord = 2; % y (up-down)
feat = 4; % jaw

dlc = obj.traj{view}; % cell array of dlc trajectories for each trial

%% first need to get jaw phase

% jawphase = warp_calculateJawPhase_localmin(dlc,feat,view,obj.bp.ev.goCue);
jawphase = warp_calculateJawPhase_zerocross(dlc,feat,view,obj.bp.ev.goCue,obj.bp.ev.lickL);


%%

nLicks = params.nLicks;
probe = prbnum;
dt = 1/400;

%%
% half rise time to first peak of jaw onset after go cue


% filter params
opts.f_cut = 60; % cutoff freq for butter filt
opts.f_n   = 2;  % filter order

% peak finding params
opts.minpkdist = 0.06; % number of ms around peaks to reject peaks
opts.minpkprom = 5;   % a threshold for the peak size


% 


% med = struct();
% jawStart = getJawTimes(view,feat,obj,opts,nLicks); % jaw oscillations proxy for lick times
% med = findMedianJawTimes(med,jawStart,nLicks);
% pfit = trialWarpFits_Jaw(jawStart,med,obj,nLicks);

% need to keep 'med' for alignSpikes
obj.med = med;


% warp spike times for each cluster (only spike times between go cue and first params.nLicks licks get
% warped)
obj = warpSpikes(obj,probe,jawStart,med,pfit);

% obj.clu{probe}(1)

% warp video data
obj = warpVideo(obj,jawStart,med,pfit);

% warp bpod lick contact times
%obj = warpLickRaster(obj,jawStart,med,pfit);  % not working



end
