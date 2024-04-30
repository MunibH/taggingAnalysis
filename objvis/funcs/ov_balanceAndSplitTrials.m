function trialsBalanced = ov_balanceAndSplitTrials(trialid, cond2use, balancegroup)
% Input:
% - trialid: Cell array where each entry is the trial IDs for a given condition
% - cond2use: Conditions to draw trials from trialid
% - balanceGroup: vector size of cond2use. conditions with same
%                 balancegroup will have their trials balance

% allTrials = cell2mat(trialid(cond2use)');

rng(1)

% Count the number of trials for each condition
trialsPerCond = trialid(cond2use);
numTrialsPerCond = cell2mat(cellfun(@numel , trialsPerCond,'uni',0));

% Get unique groups in the balancegroup vector
uniqueGroups = unique(balancegroup,'stable');


mins = [];
% Iterate through each group
for group = uniqueGroups
    % Find the indices of the cell array entries that belong to the current group
    groupIndices = balancegroup == group;

    % Get the lengths of the entries in the current group
    lengths = cellfun(@numel, trialsPerCond(groupIndices));

    % Find the minimum length in the current group
    minLength = min(lengths);

    mins = [mins ; repmat(minLength,sum(groupIndices),1)];
end

for i = 1:numel(cond2use)
    trialsBalanced{i} = trialsPerCond{i}(randperm(mins(i)));
end


end