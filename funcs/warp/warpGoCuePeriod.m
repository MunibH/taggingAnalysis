function [obj,med] = warpGoCuePeriod(obj,probe,med)
% duration of goCue to trial end period is variable
% warp this duration to median duration across all trials

% get trial end times from video data
trialEnd = getTrialEnd(obj);
med.trialEnd = median(trialEnd); % just a guess, can you video data to get actual trial end time
med.goCuePeriod = med.trialEnd - med.goCue;

for cluix = 1:numel(obj.clu{probe})
    for trix = 1:obj.bp.Ntrials
        
        x = [med.goCue trialEnd(trix)];
        y = [med.goCue med.trialEnd];
        p = polyfit(x,y,1); % linear time warp
        
        % find spikes between sample and delay to warp
        % find spike times for current trial
        spkmask = ismember(obj.clu{probe}(cluix).trial,trix);
        spkix = find(spkmask);
        spktm = obj.clu{probe}(cluix).trialtm(spkmask);
        
        mask = (spktm>=med.goCue) & (spktm<=trialEnd(trix));
        tm = spktm(mask);
        warptm = polyval(p,tm);
        obj.clu{probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
        
        % warp video data
        ft = obj.traj{1}(trix).frameTimes_warped;
        frameix = (ft>=med.goCue) & (ft<=trialEnd(trix));
        warptm = polyval(p,ft(frameix));
        obj.traj{1}(trix).frameTimes_warped(frameix) = warptm;
        
    end
end



end % warpGoCuePeriod

%%
function trialEnd = getTrialEnd(obj)
vid = obj.traj{1};
trialEnd = zeros(numel(vid),1);
for trix = 1:numel(vid)
    trialEnd(trix) = vid(trix).frameTimes(end);
end
end % getTrialEnd




















