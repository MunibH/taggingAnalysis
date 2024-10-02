function meta = tjM1_SessionMeta(meta,datapth)

meta = allSessionMeta(meta,datapth);

useSess = false(numel(meta),1);
for i = 1:numel(meta)
    r = meta(i).region;
    nProbes = numel(r);
    use = false(nProbes,1);
    for j = 1:nProbes
        if ~contains(r{j},'tjM1')
            use(j) = false;
        else
            use(j) = true;
        end
        if any(use)
            useSess(i) = true;
        else
            useSess(i) = false;
        end
    end
    meta(i).probe = meta(i).probe(use);
    meta(i).region = meta(i).region(use);
end

meta = meta(useSess);

end