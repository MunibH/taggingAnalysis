function [mask_,allquality] = findClustersByQuality(clu, metrics, desiredQuality)
 
% returns index of obj.clu{prbnum} where unit is of desiredQuality
% the desiredQuality is either based on bombcell metrics of self-labeled
% metrics from Phy. I am using bombcell metrics

%%

% assign bombcell unittype to each unit using phy_clusterID
% tagged units phy_clusterID == NaN, but gets converted to 0 for some
% reason. 
cluid = nan(numel(clu),1);
for i = 1:numel(clu)
    cluid(i) = clu(i).phy_clusterID;
end
nanmask = isnan(cluid);
cluid = cluid(~nanmask);

phyid = [metrics(:).phy_clusterID]'; 

% check if cluid and phyid are equal
if ~isequal(cluid,phyid)
    error('cluster phy ids and metrics phy ids do not match, something might be wrong with creating data object');
end

% get qualities
bcquality = {metrics(:).UnitType}';

mask_ = ismember(bcquality,desiredQuality);

% add back tagged units
mask_(end:end+sum(nanmask)) = true;

allquality = cat(1,bcquality,repmat({'tagged'},sum(nanmask),1));


%%


% % find idx where qualityList contains at least one of the patterns in
% % qualities
% 
% % remove leading/trailing spaces
% % also check if each entry in qualityList is a string/char, if not, change
% % it
% for i = 1:numel(qualityList)
%     if ~ischar(qualityList{i}) && ~isstring(qualityList{i})
%         qualityList{i} = '';
%     end
% end
% qualityList = strtrim(qualityList);
% 
% 
% if any(strcmpi(qualities,'all')) || strcmpi(qualities{1},'all')
%     idx = find(~ismember(qualityList,{'garbage','noise','trash'}));
% %     idx = 1:numel(qualityList);
%     return
% end
% 
% % handle unlabeled cluster qualities
% for i = 1:numel(qualityList)
%     if isempty(qualityList{i})
%         qualityList(i) = {'nan'};
%     end
% %     qualityList(i) = lower(qualityList(i));
% end
% 
% [~,mask] = patternMatchCellArray(qualityList, qualities, 'any');
% 
% idx = find(mask);




end % findClusters