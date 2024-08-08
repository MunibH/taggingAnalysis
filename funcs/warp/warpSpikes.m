function obj = warpSpikes(obj,probe,jawStart,med,pfit)

dt = 1/400;

for cluix = 1:numel(obj.clu{probe})
%     disp(['Warping spikes for cluster ' num2str(cluix) ' / ' num2str(numel(obj.clu{probe}))])
    obj.clu{probe}(cluix).trialtm_warped = obj.clu{probe}(cluix).trialtm;
    for trix = 1:obj.bp.Ntrials
        % find spike times for current trial
        spkmask = ismember(obj.clu{probe}(cluix).trial,trix);
        spkix = find(spkmask);
        spktm = obj.clu{probe}(cluix).trialtm(spkmask);
        
        gc = obj.bp.ev.goCue(trix);
        ls = jawStart{trix}; % current trial lick(lix) start times
        % if no licks found in current trial, move on to the next trial
        if numel(ls) == 0
            continue
        end
        
        for lix = 1:numel(ls) % lick index for current trial
            p_cur = pfit{trix,lix}; % current fit parameters for current trial and lick number

            % find spike ix b/w to warp
            if lix == 1 && lix == numel(ls) % only 1 lick
                mask = (spktm>=gc) & (spktm<=(ls(1)+median(diff(med.lickStart))));
            elseif lix == 1
                mask = (spktm>=gc) & (spktm<=(ls(2)-dt));
            elseif lix == numel(ls)
                mask = (spktm>=ls(lix)) & (spktm<=(ls(lix)+median(diff(med.lickStart))) );
            else
                mask = (spktm>=ls(lix)) & (spktm<=(ls(lix+1)-dt) );
            end
            
            tm = spktm(mask);
            
            % warp
            warptm = polyval(p_cur,tm);
            
            % assign new spike times
            obj.clu{probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
            
        end
    end
end

end % warpSpikes
