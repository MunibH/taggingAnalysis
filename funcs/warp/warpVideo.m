function obj = warpVideo(obj,jawStart,med,pfit)

dt = 1/400;

for view = 1:numel(obj.traj)

    for trix = 1:obj.bp.Ntrials
        obj.traj{view}(trix).frameTimes_warped = obj.traj{view}(trix).frameTimes - 0.5;
        
        ft = obj.traj{view}(trix).frameTimes_warped;
        gc = obj.bp.ev.goCue(trix);
        ls = jawStart{trix}; % current trial lick(lix) start times
        % if no licks found in current trial, move on to the next trial
        if numel(ls) == 0
            continue
        end
        
        for lix = 1:numel(ls) % lick index for current trial
            p_cur = pfit{trix,lix}; % current fit parameters for current trial and lick number
            
            % find frametimes to warp
            if lix == 1 && lix == numel(ls) % only 1 lick
                mask = (ft>=gc) & (ft<=(ls(1)+median(diff(med.lickStart))));
            elseif lix == 1
                mask = (ft>=gc) & (ft<=(ls(2)-dt)); %indices in frameTimes that we should warp if it's the first lick
            elseif lix == numel(ls)
                mask = (ft>=ls(lix)) & (ft<=(ls(lix)+median(diff(med.lickStart))) );
            else
                mask = (ft>=ls(lix)) & (ft<=(ls(lix+1)-dt) );
            end
            
            % warp video data
            warptm = polyval(p_cur,ft(mask));
            obj.traj{view}(trix).frameTimes_warped(mask) = warptm;
            
        end
    end
end

end % warpVideo
