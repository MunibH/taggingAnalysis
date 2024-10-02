function [mu, sd] = baselineFR(obj,params,prbnum)

% baselineFR - (nCells,1) mean FR for each clu during ITI
%%

dt = params.dt;
edges = -0.25:dt:mode(obj.bp.ev.sample);
duration = edges(end) - edges(1);

for i = 1:numel(params.cluid)
    curClu = params.cluid(i);
    trialtm = obj.clu{prbnum}(curClu).trialtm;
    trial = obj.clu{prbnum}(curClu).trial;
    for j = 1:obj.bp.Ntrials
        spkix = ismember(trial,j);        
        baselinedat(:,j,i) = histc(trialtm(spkix), edges) ./ dt; % (time,trial,unit)
    end
end

mu = squeeze(mean(baselinedat,[1,2]));
sd = squeeze(std(baselinedat,0,[1,2]));

end
