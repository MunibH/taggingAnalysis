function [obj,cluid] = QMRejection(obj,cluid,qm,prbnum)
%%

% the first two if statements here do the same thing, just separate checks
% if metrics exist and if it's empty
if ~isfield(obj,'metrics')
    % just remove based on firing rate
    meanfr = squeeze(mean(mean(obj.trialdat,1),3));
    use = ~(meanfr<qm.firing_rate);

    % always keep tagged units
    nTagged = numel(obj.tag);    
    use(end-nTagged+1:end) = true; % keep tagged units always!!

    
    cluid = cluid(use);
    obj.psth = obj.psth(:,use,:);
    obj.trialdat = obj.trialdat(:,use,:);
    return
end
if isempty(obj.metrics)
    % just remove based on firing rate
    meanfr = squeeze(mean(mean(obj.trialdat,1),3));
    use = ~(meanfr<qm.firing_rate);

    % always keep tagged units
    nTagged = numel(obj.tag);    
    use(end-nTagged+1:end) = true; % keep tagged units always!!

    
    cluid = cluid(use);
    obj.psth = obj.psth(:,use,:);
    obj.trialdat = obj.trialdat(:,use,:);
    return
end

frmask = obj.metrics{prbnum}.firing_rate > qm.firing_rate;
isimask = obj.metrics{prbnum}.isi_viol < qm.isi_viol;
presencemask = obj.metrics{prbnum}.presence_ratio > qm.presence_ratio;

use = frmask & isimask & presencemask;

% since tagged units aren't in obj.metrics, we have to add them manually
nTagged1 = numel(cluid) - numel(use);
nTagged2 = numel(obj.tag);
assert(nTagged1==nTagged2,'Number of tagged units in obj.clu and obj.tag do not match')

use(end+1:end+nTagged1) = true; % keep tagged units always!!

% remove units
cluid = cluid(use);
obj.psth = obj.psth(:,use,:);
obj.trialdat = obj.trialdat(:,use,:);


end % removeLowFRClusters