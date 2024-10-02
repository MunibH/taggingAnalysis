function [mu, sd] = trialFR(obj,params,prbnum)

% baselineFR - (nCells,1) mean FR for each clu over entire trial

dt = params.dt;
edges = 0:dt:mode(obj.bp.ev.reward) + 3;
duration = edges(end) - edges(1);

for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    trialtm = obj.clu{prbnum}(curClu).trialtm;
    trial = obj.clu{prbnum}(curClu).trial;
    for j = 1:obj.bp.Ntrials
        spkix = ismember(trial,j);        
        dat(:,j,i) = histc(trialtm(spkix), edges) ./ dt; % (time,trial,unit)
    end
end

mu = squeeze(mean(dat,[1,2]));
sd = squeeze(std(dat,0,[1,2]));

end
