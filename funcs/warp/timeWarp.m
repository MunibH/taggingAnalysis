function timeWarp(~,~,fig)
p = guidata(fig);
h = guidata(p.parentfig);

% PARSE GUI DATA
probe = h.probeList.Value;

nLicks = str2double(p.nLicks.String);

%% time warp using jaw

dt = 1/400;

% half rise time to first peak of jaw onset after go cue
view = 1; % side cam
feat = 2; % jaw

% filter params
opts.f_cut = 60; % cutoff freq for butter filt
opts.f_n   = 2;  % filter order

% peak finding params
opts.minpkdist = 0.06; % number of ms around peaks to reject peaks
opts.minpkprom = 10;   % a threshold for the peak size

med = struct();
% % warp sample and delay periods to median times (only useful if want to
% % look at early trials as defined by bpod)
% [h.obj,med] = warpSamplePeriod(h.obj,probe);
% [h.obj,med] = warpDelayPeriod(h.obj,probe,med);
% [h.obj,med] = warpGoCuePeriod(h.obj,probe,med);

jawStart = getJawTimes(view,feat,h.obj,opts,nLicks); % jaw oscillations proxy for lick times
med = findMedianJawTimes(med,jawStart,nLicks);
pfit = trialWarpFits_Jaw(jawStart,med,h.obj,nLicks);

% warp spike times for each cluster (only spike times between go cue and first params.nLicks licks get
% warped)
h.obj = warpSpikes(h.obj,probe,jawStart,med,pfit);

% warp video data
h.obj = warpVideo(h.obj,jawStart,med,pfit);

% warp bpod lick contact times
h.obj = warpLickRaster(h.obj,jawStart,med,pfit);

%% old time warping

% % get first params.Nlicks lick times, for each trial
% [lickStart,lickEnd,lickDur] = findLickTimes(h.obj,nLicks);
% 
% %----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% %find median lick time for each lick across trials
% % only using trials where a lick was present in median calculation
% med = findMedianLickTimes(lickStart,lickEnd,lickDur, nLicks);
% 
% %----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% % find fit for each trial and each lick
% pfit = trialWarpFits(lickStart,lickEnd,lickDur,med,h.obj,nLicks);
% 
% %----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----%----
% % warp spike times for each cluster (only spike times between go cue and first params.nLicks licks get
% % warped)
% for cluix = 1:numel(h.obj.clu{probe})
%     
%     disp(['Warping spikes for cluster ' num2str(cluix) ' / ' num2str(numel(h.obj.clu{probe}))])
%     h.obj.clu{probe}(cluix).trialtm_warped = h.obj.clu{probe}(cluix).trialtm;
% 
%     for trix = 1:h.obj.bp.Ntrials
%         % find spike times for current trial
%         spkmask = ismember(h.obj.clu{probe}(cluix).trial,trix);
%         spkix = find(spkmask);
%         spktm = h.obj.clu{probe}(cluix).trialtm(spkmask);
%         
%         gc = h.obj.bp.ev.goCue(trix);
%         ls = lickStart{trix}; % current trial lick(lix) start times
%         le = lickEnd{trix}; % current trial lick(lix) end times
%         ld = lickDur{trix}; % current trial lick(lix) durations
%                 
%         for lix = 1:numel(lickStart{trix}) % lick index for current trial
%             p_cur = pfit{trix,lix}; % current fit parameters for current trial and lick number
%             
%             % find spike ix b/w to warp
%             if lix==1
%                 mask = (spktm>=gc) & (spktm<=le(lix));
%             else
%                 mask = (spktm>le(lix-1)) & (spktm<=le(lix));
%             end
%             
%             tm = spktm(mask);
%             
%             % warp
%             warptm = polyval(p_cur,tm);
%             h.obj.clu{probe}(cluix).trialtm_warped(spkix(mask)) = warptm;
%         end
%         
%     end
%     
% end

%% mark that we've warped data, set warp checkbox to visible and checked,
% turn off others
h.warped = 1;
h.psthDataList.Value = 3;

f = msgbox('Warping Completed');
pause(1);
delete(f);

close(figure(535))

guidata(h.fig(1), h);

probe = get(h.probeList, 'Value');
unit = get(h.unitList, 'Value');

clu = h.obj.clu{probe}(unit);

updateRaster(h.fig(1), clu);
updatePSTH(h.fig(1), clu);
updateBehav([],[],h.fig(1));
updateVideo([],[],h.fig(1))


end % timeWarp
















