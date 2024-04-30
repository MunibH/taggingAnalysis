function quality = getUnitQuality(clu,cluid)
% get qualities from clu, for cluster ids in cluid

quality = {clu(cluid).quality};
for i = 1:numel(quality)
    if ~ischar(quality{i}) && ~isstring(quality{i})
        quality{i} = '';
    end
end
quality = strtrim(quality);

end