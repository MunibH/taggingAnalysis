function params = defaultParams()

params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0.1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for (and get trial numbers for)
params.condition(1)         = {'(hit|miss|no)'};                             % all trials
params.condition(end+1)     = {'R&hit&~autowater'};                                     % right hit trials
params.condition(end+1)     = {'L&hit&~autowater'};                                     % left hit trials
params.condition(end+1)     = {'R&miss&~autowater'};                                     % right hit trials
params.condition(end+1)     = {'L&miss&~autowater'};                                     % left hit trials

% time from align event to grab data for
params.tmin = -2.1;
params.tmax = 2;
params.dt = 1/100; % bin size

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect';

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.advance_movement = 0.0; % not sure if still being used anywhere, leaving in place just in case

params.behav_only = 1; % if 1, don't process neural data - speeds up analysis if you only need behavior data

params.remove_tag = 1; % if 1, removes tagging trials before doing anything else

params.region = {'ALM','tjM1','MC'}; % regions to use

params.events = {'bitStart','sample','delay','goCue'};

params.qm.presence_ratio = 0.9; % greater than
params.qm.firing_rate = 1; % greater than
params.qm.isi_viol = 10; % less than
%'unit_amp > 150'

end % defaultParams()