function obj = TrimDataEarlyLickTrialsSampleANDDelay(obj,samp_start,samp_end,delay_start,delay_end,trial)


%% ME and TRAJ

traj1 = obj.traj{1}(trial);
traj2 = obj.traj{2}(trial);
me = obj.me{trial};
ft = traj1.frameTimes - 0.5;

[ft_adjusted,ft_mask,timeAdjustments] = EarlyLickFormatFrameTimes(ft, {[samp_start,samp_end];[delay_start,delay_end]});

% figure; hold on;
% plot(ft,obj.traj{1}(trial).ts(:,[1,2],4))
% xline(obj.bp.ev.goCue(trial))

obj.traj{1}(trial).frameTimes = ft_adjusted;
obj.traj{1}(trial).ts = obj.traj{1}(trial).ts(ft_mask,:,:);
obj.me{trial} = obj.me{trial}(ft_mask);

% figure; hold on;
% plot(obj.traj{1}(trial).frameTimes, obj.traj{1}(trial).ts(:,[1,2],4))
% xline(mode(obj.bp.ev.goCue))

%% EVENTS

obj.bp.ev.delay(trial) = mode(obj.bp.ev.delay);
obj.bp.ev.goCue(trial) = mode(obj.bp.ev.goCue);
obj.bp.ev.reward(trial) = obj.bp.ev.reward(trial) - sum(timeAdjustments);

tempL = obj.bp.ev.lickL{trial};
if ~isempty(tempL)
    mask = obj.bp.ev.lickL{trial}>=samp_start;
    tempL(mask) = tempL(mask) - timeAdjustments(1);
    mask = obj.bp.ev.lickL{trial}>=delay_start;
    tempL(mask) = tempL(mask) - timeAdjustments(2);
    obj.bp.ev.lickL{trial} = tempL;
end


tempR = obj.bp.ev.lickR{trial};
if ~isempty(tempR)
    mask = obj.bp.ev.lickR{trial}>=samp_start;
    tempR(mask) = tempR(mask) - timeAdjustments(1);
    mask = obj.bp.ev.lickR{trial}>=delay_start;
    tempR(mask) = tempR(mask) - timeAdjustments(2);
    obj.bp.ev.lickR{trial} = tempR;
end

%% CLU

for iprobe = 1:numel(obj.clu)
    clu = obj.clu{iprobe};
    for iclu = 1:numel(clu)
        c = clu(iclu);
        c_adjusted = EarlyLickFormatClu(c,{[samp_start,samp_end];[delay_start,delay_end]},trial);
        obj.clu{iprobe}(iclu) = c_adjusted;
    end
end

end