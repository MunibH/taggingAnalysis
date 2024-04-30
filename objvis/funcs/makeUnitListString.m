function unitListString = makeUnitListString(h)

nUnits = numel(h.params.cluid);

% get unit qualities 
quality = getUnitQuality(h.obj.clu{1},h.params.cluid);

% get num spikes
for i = 1:nUnits
    nspks(i) = numel(h.obj.clu{1}(h.params.cluid(i)).tm);
end

% make unit string
unitListString = arrayfun(@(i) sprintf('Unit %d, %s spks, %s', i, num2str(nspks(i)), quality{i}), 1:numel(h.params.cluid), 'UniformOutput', false);

end