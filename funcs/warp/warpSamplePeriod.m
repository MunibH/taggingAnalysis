function [obj,med] = warpSamplePeriod(obj,probe)
% start of sample occurs at same time in trial
% duration of sample to delay period is variable
% warp this duration to median duration across all trials
med.sample = median(obj.bp.ev.sample);
med.delay = median(obj.bp.ev.delay);
med.samplePeriod = med.delay - med.sample;
obj.bp.ev.delay_warped = obj.bp.ev.delay;

for cluix = 1:numel(obj.clu{probe})
    obj.clu{probe}(cluix).trialtm_warped = obj.clu{probe}(cluix).trialtm;
    for trix = 1:obj.bp.Ntrials
        obj.traj{1}(trix).frameTimes = obj.traj{1}(trix).frameTimes;
        obj.traj{1}(trix).frameTimes_warped = obj.traj{1}(trix).frameTimes - 0.5;  % correct for 0.5 sec shift
        if obj.bp.ev.delay(trix) ~= med.delay
            delay = obj.bp.ev.delay(trix); % when does delay period start
            obj.bp.ev.delay_warped(trix) = med.delay;
            x = [med.sample delay];
            y = [med.sample med.delay];
            p = polyfit(x,y,1); % linear time warp
            
            % find spikes between sample and delay to warp
            % find spike times for current trial
            spkmask = ismember(obj.clu{probe}(cluix).trial,trix);
            spkix = find(spkmask);
            spktm = obj.clu{probe}(cluix).trialtm(spkmask);
            
            mask = (spktm>=med.sample) & (spktm<=delay);
            tm = spktm(mask);
            warptm = polyval(p,tm);
            obj.clu{probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
            
            % warp video data
            ft = obj.traj{1}(trix).frameTimes_warped;
            frameix = (ft>=med.sample) & (ft<=delay);
            warptm = polyval(p,ft(frameix));
            obj.traj{1}(trix).frameTimes_warped(frameix) = warptm;
        end
    end
end



end % warpSamplePeriod






