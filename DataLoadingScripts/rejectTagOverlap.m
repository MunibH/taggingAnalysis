function use = rejectTagOverlap(obj,thresh,prbnum)

% returns logical mask size of obj.clu{prbnum}. 1 means use that unit in
% analysis, 0 don't use

use = false(numel(obj.clu{prbnum}),1); % all units, including tagged

% get tagged units for this probe
tagprobe = [obj.tag(:).probenum]';
usetag = tagprobe==prbnum;
tag = obj.tag(usetag);

nTag = numel(tag);
nClu = numel(obj.clu{prbnum});
nClu_min_nTag = nClu - nTag;
use(nClu_min_nTag+1:end) = true; % keep tagged units always!!


for itag = 1:numel(tag)
    overlap(:,itag) = tag(itag).tagOverlap(1:nClu_min_nTag) > thresh;
end

discard = sum(overlap,2) > 0;
use(1:nClu_min_nTag) = ~discard;




end % rejectTagOverlap