function [lick] = GetLickPhase(obj,p)

% Normalized cutoff frequencies (between 0 and 1)
Wn = [p.filt.Fc1 p.filt.Fc2] / (p.filt.Fs / 2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(4, Wn, 'bandpass');


traj = obj.traj{1}; % side cam
feat = 'jaw';
ifeat = ismember(traj(1).featNames,feat);

lickct = 1;
for itrial = 1:obj.bp.Ntrials
    % get lick times for current trial
    licktimes = sort([obj.bp.ev.lickL{itrial} obj.bp.ev.lickR{itrial}]);

    % get jaw position for current trial
    jaw = traj(itrial).ts(:,2,ifeat); % jaw, side cam, ycoord
    % frame times
    ft = traj(itrial).frameTimes - 0.5;
    
    % only keep lick times that are on video
    licktimes(licktimes<ft(1) | licktimes>ft(end)) = [];

    % loop through each bpod lick
    for ilick = 1:numel(licktimes)

        % find index in frame times where this lick occurred
        tix = findTimeIX(ft,licktimes(ilick));
        % get surrounding presamp and postsamp indicies
        startix = max(1,tix-p.presamp);
        endix = min(tix+p.postsamp,numel(ft));

        % extract the lick
        x = jaw(startix:endix);

        % find start and end time of 'lick'
        licktimes_ft = ft(startix:endix);


        % bandpass filter the lick (fill missing data first)
        x_filled = fillmissing(x, 'movmean', p.fillmiss_window);  % Replace window_size with an appropriate value
        if any(isnan(x_filled)) % there was one lick I found that fillmissing still didn't fill all values, skip that one...
            continue
        end
        lick_filt = filtfilt(b,a,x_filled);

        % where bpod detected lick
        center_idx = round(numel(lick_filt)/2);
        % correct the center idx (bpod lick time might not actually be peak of lick
        [~,center_idx_] = max(lick_filt(center_idx-p.peri_center:center_idx+p.peri_center));
        center_idx_corrected = center_idx_ + center_idx - p.peri_center;
        % walk backward from center find start of lick
        [~, min_back_idx] = min(lick_filt(1:center_idx_corrected));  % Index of minimum before the peak
        % walk forward from center find end of lick
        [~, min_forward_idx] = min(lick_filt(center_idx_corrected:end));  % Index of minimum after the peak
        min_forward_idx = min_forward_idx + center_idx_corrected - 1;  % Adjust index relative to the original signal

        % store lick and velocity
        lick.pos{lickct} = lick_filt(min_back_idx:min_forward_idx);
        lick.vel{lickct} = gradient(lick.pos{lickct});

        % don't think I'm suing this anymore
        % ix.lick_center{lickct} = center_idx_corrected;


        % now I have the lick, I need the lick start and end times within
        % the trial (frameTimes -- ft)
        lick.frame_times{lickct} = licktimes_ft(min_back_idx:min_forward_idx);
        lick.trial(lickct) = itrial;
        lick.trialtm(lickct) = licktimes_ft(center_idx_corrected);
        lickct = lickct + 1;

    end


end
lick.trial = lick.trial';
lick.pos = lick.pos';
lick.vel = lick.vel';
lick.frame_times = lick.frame_times';
lick.trialtm = lick.trialtm';
lick.nLicks = numel(lick.pos);


% get lick time series, padded out to length of max num time points in lick.pos
% get lick phase time series at same time
durations = cell2mat(cellfun(@numel,lick.pos,'uni',0));
maxLickDuration = max(durations);
medianLickDuration = median(durations);
lick.pos_ts = nan(maxLickDuration,lick.nLicks);
lick.vel_ts = nan(maxLickDuration,lick.nLicks);
lick.pos_phase = nan(medianLickDuration,lick.nLicks);
lick.vel_phase = nan(medianLickDuration,lick.nLicks);
for iLick = 1:lick.nLicks
    d = durations(iLick);
    lick.pos_ts(1:d,iLick) = lick.pos{iLick};
    lick.vel_ts(1:d,iLick) = lick.vel{iLick};

    % phase
    % for each lick, assign linspace(0,2*pi,duration)
    thisphase = linspace(0,2*pi,d);
    % find ix where lick trialtm falls in frame_times, assign lick.phase as thisphase(ix)
    ix = findTimeIX(lick.frame_times{iLick},lick.trialtm(iLick));
    lick.peak_pro_phase(iLick) = thisphase(ix);

    % for time series, downsample to shortest lick for each of the licks
    lick.pos_phase(:,iLick) = resample(lick.pos{iLick}, medianLickDuration, d);
    lick.vel_phase(:,iLick) = resample(lick.vel{iLick}, medianLickDuration, d);

end
lick.peak_pro_phase = lick.peak_pro_phase';
lick.phase_ts_x = linspace(0, 2*pi, medianLickDuration)';  % Assign phases from 0 to 2*pi


end