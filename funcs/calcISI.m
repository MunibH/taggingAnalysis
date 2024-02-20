function [nISI, ISIcoord, mx] = calcISI(obj,probenum,ISIedges)

allclu = obj.clu{probenum};

for i = 1:numel(allclu)
    clu = allclu(i);
    
    dt = diff(clu.tm);
    dt = dt(diff(clu.trial)==0);
    if isempty(dt)
        dt = 1;
    end
    
    nISI(:,i) = histc([dt; -dt], ISIedges);
    ISIcoord(:,i) = ISIedges+mean(diff(ISIedges))./2;
    
    mx(i) = max(nISI(:,i));
    mx(i) = 1.05.*mx(i);
    
end

end