function use = QMRejection(obj,qm,prbnum)

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

% firing rate
% get number of spikes per cluster
for i = 1:nClu_min_nTag
    nspks(i) = numel(obj.clu{prbnum}(i).tm);
end
if isfield(obj.sglx,'imec')
    end_time = (obj.sglx(prbnum).imec.fileStart(end)) + (obj.sglx(prbnum).imec.Nsamp(end)./obj.sglx(prbnum).imec.fs);
    duration =  end_time - obj.sglx(prbnum).imec.fileStart(1); % in sec
else
    end_time = (obj.sglx(prbnum).fileStart(end)) + (obj.sglx(prbnum).Nsamp(end)./obj.sglx(prbnum).fs);
    duration =  end_time - obj.sglx(prbnum).fileStart(1); % in sec
end
firing_rate = nspks./duration;
frmask = firing_rate >= qm.firing_rate;
frmask(end+1:numel(use)) = true; % keep tagged units always!!

isimask = [obj.metrics{prbnum}(:).fractionRPVs_estimatedTauR] < qm.isi_viol;
isimask(end+1:numel(use)) = true; % keep tagged units always!!

presencemask = [obj.metrics{prbnum}(:).presenceRatio] > qm.presence_ratio;
presencemask(end+1:numel(use)) = true; % keep tagged units always!!

use = (frmask & isimask & presencemask)';



%
% % the first two if statements here do the same thing, just separate checks
% % if metrics exist and if it's empty
% if ~isfield(obj,'metrics')
%     % just remove based on firing rate
%     meanfr = squeeze(mean(mean(obj.trialdat,1),3));
%     use = ~(meanfr<qm.firing_rate);
%
%     % always keep tagged units
%     nTagged = numel(obj.tag);
%     use(end-nTagged+1:end) = true; % keep tagged units always!!
%
%
%     cluid = cluid(use);
%     obj.psth = obj.psth(:,use,:);
%     obj.trialdat = obj.trialdat(:,use,:);
%     return
% end
% if isempty(obj.metrics)
%     % just remove based on firing rate
%     meanfr = squeeze(mean(mean(obj.trialdat,1),3));
%     use = ~(meanfr<qm.firing_rate);
%
%     % always keep tagged units
%     nTagged = numel(obj.tag);
%     use(end-nTagged+1:end) = true; % keep tagged units always!!
%
%
%     cluid = cluid(use);
%     obj.psth = obj.psth(:,use,:);
%     obj.trialdat = obj.trialdat(:,use,:);
%     return
% end
% %%
%
% if numel(fieldnames(obj.metrics{prbnum})) > 5 % bombcell
%     use = false(numel(obj.clu{prbnum}),1); % all units, including tagged
%
%     nTag = numel(obj.tag);
%     nClu = numel(obj.clu{prbnum});
%     nClu_min_nTag = nClu - nTag;
%     warning('QMRejection only works for 1 probe recordings currently')
%     warning('It assumes that all tagged units come from all probes')
%     use(nClu_min_nTag+1:end) = true; % keep tagged units always!!
%
%     % firing rate
%     % get number of spikes per cluster
%     for i = 1:nClu_min_nTag
%         nspks(i) = numel(obj.clu{prbnum}(i).tm);
%     end
%     end_time = (obj.sglx(prbnum).imec.fileStart(end)) + (obj.sglx(prbnum).imec.Nsamp(end)./obj.sglx(prbnum).imec.fs);
%     duration =  end_time - obj.sglx(prbnum).imec.fileStart(1); % in sec
%     firing_rate = nspks./duration;
%     frmask = firing_rate >= qm.firing_rate;
%     frmask(end+1:numel(use)) = true; % keep tagged units always!!
%
%     isimask = [obj.metrics{prbnum}(:).fractionRPVs_estimatedTauR] < qm.isi_viol;
%     isimask(end+1:numel(use)) = true; % keep tagged units always!!
%
%     presencemask = [obj.metrics{prbnum}(:).presenceRatio] > qm.presence_ratio;
%     presencemask(end+1:numel(use)) = true; % keep tagged units always!!
%
%     use = frmask & isimask & presencemask;
%
% else % custom quality metrics before using bombcell
%     frmask = obj.metrics{prbnum}.firing_rate >= qm.firing_rate;
%     isimask = obj.metrics{prbnum}.isi_viol < qm.isi_viol;
%     presencemask = obj.metrics{prbnum}.presence_ratio > qm.presence_ratio;
%
%     use = frmask & isimask & presencemask;
%
%     % since tagged units aren't in obj.metrics, we have to add them manually
%     nTagged1 = numel(cluid) - numel(use);
%     nTagged2 = numel(obj.tag);
%     assert(nTagged1==nTagged2,'Number of tagged units in obj.clu and obj.tag do not match')
%
%     use(end+1:end+nTagged1) = true; % keep tagged units always!!
% end


% %%
%
% % remove units
% cluid = cluid(use);
% obj.psth = obj.psth(:,use,:);
% obj.trialdat = obj.trialdat(:,use,:);


end % QMRejection