function [obj,med] = warpDelayPeriod(obj,probe,med)
% duration of delay to goCue period is variable
% warp this duration to median duration across all trials
med.goCue = median(obj.bp.ev.goCue);
med.delayPeriod = med.goCue - med.delay;
obj.bp.ev.goCue_warped = obj.bp.ev.goCue;

for cluix = 1:numel(obj.clu{probe})
    for trix = 1:obj.bp.Ntrials
        if obj.bp.ev.goCue(trix) ~= med.goCue
            goCue = obj.bp.ev.goCue(trix); % when does goCue period start
            obj.bp.ev.goCue_warped(trix) = med.goCue;
            x = [med.delay goCue];
            y = [med.delay med.goCue];
            p = polyfit(x,y,1); % linear time warp
            
            % find spikes between sample and delay to warp
            % find spike times for current trial
            spkmask = ismember(obj.clu{probe}(cluix).trial,trix);
            spkix = find(spkmask);
            spktm = obj.clu{probe}(cluix).trialtm(spkmask);
            
            mask = (spktm>=med.delay) & (spktm<=goCue);
            tm = spktm(mask);
            warptm = polyval(p,tm);
            obj.clu{probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
            
            % warp video data
            ft = obj.traj{1}(trix).frameTimes_warped;
            frameix = (ft>=med.delay) & (ft<=goCue);
            warptm = polyval(p,ft(frameix));
            obj.traj{1}(trix).frameTimes_warped(frameix) = warptm;
            
        end
    end
end



end % warpDelayPeriod






