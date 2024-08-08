function p = trialWarpFits_Jaw(jawStart,med,obj,nLicks)

dt = 1/400;

p = cell(obj.bp.Ntrials,nLicks);
mls = med.lickStart;

for trix = 1:obj.bp.Ntrials
    ls = jawStart{trix}; % current trial lick start times
        
    for lix = 1:numel(ls) % lick index
        % median values for current lick
        
        % warp lick times
        % if first lick: warp from go cue to just before 2nd lick start
        % if lick2 or lick>2: warp from start of current lick to just
        % before start of next lick
        % if last lick: warp from start of last lick to median(diff(mls))
        % [lickRate]
        if lix == 1 && lix == numel(ls) % only 1 lick
            x = [obj.bp.ev.goCue(trix) ls(1)+median(diff(mls))];
            y = [mode(obj.bp.ev.goCue) mls(1)+median(diff(mls))];

        elseif lix == 1
            x = [obj.bp.ev.goCue(trix) ls(2)-dt];
            y = [mode(obj.bp.ev.goCue) mls(2)-dt];

        elseif lix == numel(ls)
            x = [ls(lix) ls(lix)+median(diff(mls))];
            y = [mls(lix) mls(lix)+median(diff(mls))];
        else
            x = [ls(lix) ls(lix+1)-dt]; 
            y = [mls(lix) mls(lix+1)-dt];
        end
        
        % fit original time to warped time
        p{trix,lix} = polyfit(x,y,1);
        
    end
    
end

end