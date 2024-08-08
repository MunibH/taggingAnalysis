function jawStart = getJawTimes(view,feat,obj,opts,nLicks)

jawStart = cell(obj.bp.Ntrials, 1);


%%

% figure;
for trix = 1:obj.bp.Ntrials
    tm = obj.traj{view}(trix).frameTimes;
    fs = 1./median(diff(tm));
    vidx = obj.traj{view}(trix).ts(:, 1, feat);
    vidy = obj.traj{view}(trix).ts(:, 2, feat);
    
    % get jaw pos relative to origin, fill nans, mean center
    jawPos = sqrt((vidx-nanmean(vidx)).^2 + (vidy-nanmean(vidy)).^2);
    jawPos = vidy; % HACK MNE
    jawPos = fillmissing(jawPos, 'nearest');
    jawPos = jawPos-mean(jawPos);
    
    % filter jaw position
    Wn = opts.f_cut/fs/2;
    [b, a] = butter(opts.f_n, Wn);
    filtJawPos = filtfilt(b, a, jawPos);
    
    % find peaks in jaw position with min peak dist and min pk prominance
    [~, locs] = findpeaks(filtJawPos, 'MinPeakDistance', ceil(opts.minpkdist*fs), 'MinPeakProminence',opts.minpkprom);

    jtms = tm(locs);
    
    goCue = obj.bp.ev.goCue(trix)+0.5;
    

    %firstPeakIX = find(jtms>goCue+0.05, 1, 'first');
    
    peaktm = jtms((jtms>goCue+0.05));
    
    %     figure; plot(tm,filtJawPos); hold on; plot(tm(locs),filtJawPos(locs),'r.','MarkerSize',10)
    %     plot(peaktm,20*ones(size(peaktm)),'m.','MarkerSize',15)
    
    % keep first nLicks licks after go cue
    if numel(peaktm) > nLicks
        peaktm = peaktm(1:nLicks);
    end
    jawStart{trix} = peaktm - 0.5; % account for 0.5 second shift from neural to vid data (not sure why it exists)
   
    
end


end % get jawTimes