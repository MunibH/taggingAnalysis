function tag = getTagFromObj(obj,params,meta)

% get only tagging data from processed obj
for iprb = 1:numel(params.probe)
    tagix = find(ismember(params.quality{iprb},'tagged'));

    tag.nTag(iprb) = numel(tagix);
    tag.id.clu{iprb} = params.cluid{iprb}(tagix);
    tag.id.obj{iprb} = tagix;
    tag.shank{iprb} = params.shank{iprb}(tagix);
    tag.channel{iprb} = params.channel{iprb}(tagix);
    tag.region{iprb} = params.region{iprb}(tagix);
    tag.time = obj.time;
    tag.psth{iprb} = obj.psth{iprb}(:,tagix,:);
    tag.trialdat{iprb} = obj.trialdat{iprb}(:,tagix,:);
    
    tag.trialid = params.trialid;
    tag.bp = obj.bp;
    tag.anm = meta.anm;
    tag.date = meta.date;
    tag.probeType = meta.probeType;
    tag.eventTimes = params.eventTimes;

end

end