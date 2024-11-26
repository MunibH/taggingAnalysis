function facemap = LoadAndFormatFacemapPredictions(datapth,meta,obj,params)

load(fullfile(datapth, [meta.anm '_' meta.date '_FacemapPredictions.mat' ]))
% pred size = (time*trials,units)
% neuralActivity size = (time*trials,units)
% motionSVD size = (time*trials,SVs)

%% put data into trials, align to go cue


tmin = params.tmin;
tmax = params.tmax;
nBefore = abs(round(params.tmin./dt));
nAfter = round(params.tmax./dt);
nCumTime = cumsum(double(nTimeEachTrial));
nTrials = obj.bp.Ntrials;
nUnits = size(neuralActivity,2);
nSVs = size(motionSVD,2);

data.n = nan(nBefore+nAfter+1,nTrials,nUnits);
data.pred = nan(nBefore+nAfter+1,nTrials,nUnits);
data.motionSVD = nan(nBefore+nAfter+1,nTrials,nSVs);
fns = fieldnames(data);
for i = 1:nTrials

    if i == 1
        ix = 1:nCumTime(i);
    else
        ix = nCumTime(i-1)+1:nCumTime(i);
    end

    % % get go cue time
    gocue = obj.bp.ev.goCue(i);
    trialtime = (1:numel(ix)).*dt; % time in trial, aligned to go cue
    igocue = findTimeIX(trialtime,gocue);

    temp.n = neuralActivity(ix,:);
    temp.pred = pred(ix,:);
    temp.motionSVD = motionSVD(ix,:);
    
    % Determine the valid range of indices to extract
    start_idx = max(1, igocue - nBefore);
    end_idx = min(numel(ix), igocue + nAfter);

    % Determine the corresponding range in the output array
    ostart = nBefore - (igocue - start_idx) + 1;
    oend = nBefore + (end_idx - igocue) + 1;

    % data
    
    for f = 1:numel(fns)
        t = fns{f};
        data.(t)(ostart:oend,i,:) = temp.(t)(start_idx:end_idx,:);
    end

end

for f = 1:numel(fns)
    t = fns{f};
    facemap.(t) = data.(t)(1:end-1,:,:);
end


end