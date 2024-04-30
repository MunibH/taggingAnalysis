function meta = subsetMetaByParams(meta,params)
% region check
regions = lower({meta(:).region});
if strcmpi(params.region,'any')
    regioncheck = true(size(meta));
else
    regioncheck = contains(regions,lower(params.region));
end
if sum(regioncheck)==0
    error('No sessions found with specified region');
end

% probe type check
probetypes = lower({meta(:).probeType});
if strcmpi(params.probeType,'any')
    probecheck = true(size(meta));
else
    probecheck = contains(probetypes,lower(params.probeType));
end
if sum(probecheck&regioncheck)==0
    error('No sessions found with specified region and probe type');
end

meta = meta(probecheck&regioncheck);

end