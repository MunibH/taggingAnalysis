function trialNums = findTrials(obj, conditions)

% % older data objects have obj.bp.autowater.nums
% % newer data objects do not (from summer2021 pipeline)
% if ~isfield(obj.bp.autowater, 'nums')
%     tmp = obj.bp.autowater;
%     obj.bp = rmfield(obj.bp, 'autowater');
%     obj.bp.autowater.nums = tmp + (tmp-1)*-2;
% end

varnames = getStructVarNames(obj);
Ntrials = obj.bp.Ntrials;
for i = 1:numel(varnames)
    eval([varnames{i} ' = obj.bp.' varnames{i} ';']);
    
    if eval(['numel(' varnames{i} ')==obj.bp.Ntrials && isrow(' varnames{i} ')'])
        eval([varnames{i} '=' varnames{i} ''';']);
    end
    
    eval([varnames{i} '=' varnames{i} '(~isnan(' varnames{i} '))' ''';'])
end


mask = zeros(obj.bp.Ntrials, numel(conditions));

for i = 1:numel(conditions)
    try
        mask(:,i) = eval(conditions{i});
    catch
        mask(:,i) = eval(conditions{i}{1});
    end
    trialNums{i} = find(mask(:,i));
end



end % findTrials