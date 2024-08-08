function [channel,shank] = getShankAndChannels(obj,prbnum,cluid)

channel = nan(numel(cluid),1);
shank = nan(numel(cluid),1);
for iclu = 1:numel(cluid)
    chans = obj.clu{prbnum}(cluid(iclu)).channel;
    if numel(chans)>1
         channel(iclu) = mode(chans);
    else
        channel(iclu) = chans;
    end
end


if isfield(obj.clu{prbnum},'shank')
    shank = [obj.clu{prbnum}(cluid).shank]';
else % h2 probes (32 channels per shank)
    mask = (channel/32) <= 1;
    shank(mask) = 0;
    shank(~mask) = 1;
end

end