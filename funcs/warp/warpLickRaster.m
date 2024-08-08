function obj = warpLickRaster(obj,jawStart,med,pfit)

dt = 1/400;

obj.bp.ev.lickR_warped = obj.bp.ev.lickR;
obj.bp.ev.lickL_warped = obj.bp.ev.lickL;

for trix = 1:obj.bp.Ntrials
    
    gc = obj.bp.ev.goCue(trix);
    ls = jawStart{trix}; % current trial lick(lix) start times
    % if no licks found in current trial, move on to the next trial
    if numel(ls) == 0
        continue
    end
    
    for lix = 1:numel(ls) % lick index for current trial
        p_cur = pfit{trix,lix}; % current fit parameters for current trial and lick number
        
        if (obj.bp.R(trix)&&obj.bp.hit(trix)&&~obj.bp.early(trix))
            type = 'lickR';
            warptype = 'lickR_warped';
            try
                lickR = obj.bp.ev.lickR{trix}(lix);
            catch  % if more licks found than lick contacts
                break
            end
        elseif (obj.bp.hit(trix)) % HACK YH (only lickL exist in the foraging task)
            %(obj.bp.L(trix)&&obj.bp.hit(trix)&&~obj.bp.early(trix))
            type = 'lickL';
            warptype = 'lickL_warped';
            try
                lickL = obj.bp.ev.lickL{trix}(lix);
            catch % if more licks found than lick contacts
                break
            end
        else
            break
        end
        
        % find licks b/w to warp
        mask = lix;
        
        % warp video data
        warptm = polyval(p_cur,obj.bp.ev.(type){trix}(mask));
        obj.bp.ev.(warptype){trix}(mask) = warptm;        
    end
end

end % warpRaster
