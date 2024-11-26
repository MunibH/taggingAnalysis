function clu = EarlyLickFormatClu(clu,intervals,currentTrial)

intervals = sortrows(cell2mat(intervals), 1);
OnTrial = clu.trial==currentTrial;
timeAdjustments = diff(intervals');

trial = clu.trial;
trialtm = clu.trialtm;

trialtm_adjusted = trialtm;
trial_adjusted = trial;

% spikes to remove on current trial (mask(:,1) = sample   mask(:,2) = delay)
for i = 1:size(intervals, 1)
    mask(:,i) = (trialtm >= intervals(i, 1) & trialtm <= intervals(i, 2)) & OnTrial;
end

allmask = any(mask,2); % all spikes to remove across epochs

for i = 1:size(intervals,1)
    epoch_mask = (trialtm >= intervals(i,1)) & OnTrial;
    trialtm_adjusted(epoch_mask) = trialtm_adjusted(epoch_mask) - timeAdjustments(i);
end

% samp_mask = (trialtm >= intervals(1,1)) & OnTrial;
% trialtm_adjusted(samp_mask) = trialtm_adjusted(samp_mask) - timeAdjustments(1);
% delay_mask = (trialtm >= intervals(2,1)) & OnTrial;
% trialtm_adjusted(delay_mask) = trialtm_adjusted(delay_mask) - timeAdjustments(2);

trialtm_adjusted = trialtm(~allmask);
trial_adjusted = trial(~allmask);

% tempL = obj.bp.ev.lickL{trial};
% mask = obj.bp.ev.lickL{trial}>=intervals(1,1); % samp_start
% tempL(mask) = tempL(mask) - timeAdjustments(1);
% mask = obj.bp.ev.lickL{trial}>=intervals(2,1); % delay_start
% tempL(mask) = tempL(mask) - timeAdjustments(2);

% figure; plot(trialtm,trial,'.'); hold on; plot(trialtm_adjusted,trial_adjusted,'r.')

clu.trial = trial_adjusted;
clu.trialtm = trialtm_adjusted;




end
