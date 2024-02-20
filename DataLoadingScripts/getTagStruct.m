function tag = getTagStruct(meta,obj,params)

for isess = 1:numel(meta)
    % remove tagging trials
    obj(isess) = deleteTaggingTrials(obj(isess));
    % trials by condition
    tag(isess).trialid = findTrials(obj(isess), params.condition);
    params.trialid = tag(isess).trialid;

    % find tagged units
    tag(isess).cluid = [obj(isess).tag(:).cluid];
    % % tag.probeid = 1; % I don't have any sessions with dual probes, so
    % you'll have to fill this in to find the appropriate tagged unit
    tag(isess).region = meta(isess).region;

    obj(isess).clu{end+1} = obj(isess).clu{1}(tag(isess).cluid); % add a fake 'probe' of just tagged units,
    % can also just overwrite the probe if you want since this won't be saved to the data obj

    obj(isess) = alignSpikes(obj(isess),params,numel(obj(isess).clu)); % just align the new 'fake' probe data

    params.cluid = cell(numel(obj(isess).clu),1);
    params.cluid{end} = 1:numel(obj(isess).clu{end}); % process all tagged units
    temp = getSeq(obj(isess),params,numel(obj(isess).clu));

    tag(isess).time = temp.time;
    tag(isess).clu = obj(isess).clu{end};
    tag(isess).psth = temp.psth{end}; % specify probe properly
    tag(isess).trialdat = temp.trialdat{end}; % specify probe properly
    tag(isess).anm = meta(isess).anm;
    tag(isess).date = meta(isess).date;
end


end