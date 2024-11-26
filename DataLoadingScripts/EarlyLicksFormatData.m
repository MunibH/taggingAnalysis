function EarlyLicksFormatData(obj)

%%
sampleLength = mode(obj.bp.ev.delay - obj.bp.ev.sample);
delayLength = mode(obj.bp.ev.goCue - obj.bp.ev.delay);

earlyTrials = obj.bp.early;
earlyTrialsNum = find(earlyTrials);
 
% for each early lick trial, 
for i = 1:numel(earlyTrialsNum)
    itrial = earlyTrialsNum(i);
    
    % % some frametimes are messed up, skip those trials
    % traj1 = obj.traj{1}(trial);
    % ft = traj1.frameTimes - 0.5;
    % if ft(1)>0.2
    %     continue
    % end
    % 
    % if itrial == 129
    %     'a'
    % end

    % determine if early lick sample, delay, or both
    isEarlySample = round((obj.bp.ev.delay(itrial) - obj.bp.ev.sample(itrial)),4) ~= round(sampleLength,4);
    isEarlyDelay = round((obj.bp.ev.goCue(itrial) - obj.bp.ev.delay(itrial)),4) ~= round(delayLength,4);

    % get times to remove
    sampleStart = obj.bp.ev.delay(itrial) - sampleLength;
    sampleEnd = obj.bp.ev.delay(itrial);
    sampleCutStart = obj.bp.ev.sample(itrial);
    sampleCutEnd = sampleStart - 0.01;

    delayStart = obj.bp.ev.goCue(itrial) - delayLength;
    delayEnd = obj.bp.ev.goCue(itrial);
    delayCutStart = obj.bp.ev.delay(itrial);
    delayCutEnd = delayStart - 0.01;
    
    % set epoch start and end time, cut out the extra time from obj.clu,
    % obj.traj, obj.me
    if isEarlySample && isEarlyDelay
        obj = TrimDataEarlyLickTrialsSampleANDDelay(obj,sampleCutStart,sampleCutEnd,delayCutStart,delayCutEnd,itrial);
    elseif isEarlySample
        obj = TrimDataEarlyLickTrialsSampleORDelay(obj,sampleCutStart,sampleCutEnd,itrial);
    elseif isEarlyDelay
        obj = TrimDataEarlyLickTrialsSampleORDelay(obj,delayCutStart,delayCutEnd,itrial);
    end
end


obj.bp.early(:) = false;

end


