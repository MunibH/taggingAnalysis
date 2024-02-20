function [wav,chan] = getMeanWav(obj,probenum)

allclu = obj.clu{probenum};

for i = 1:numel(allclu)
    cluwav = allclu(i).spkWavs;
    
    wav(:,i) = mean(cluwav,2);
    if isfield(allclu(i),'spkWavChan')
        chan(i) = allclu(i).spkWavChan;
    elseif isfield(allclu(i),'cluster_channel')
        chan(i) = allclu(i).cluster_channel;
    else % for testing
        chan(i) = randsample(1:384,1);
    end
    
end

end