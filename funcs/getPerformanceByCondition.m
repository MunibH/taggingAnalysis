function perf = getPerformanceByCondition(obj,params)

perf = nan(numel(obj),numel(params(1).trialid)); % (sessions,conditions)
for i = 1:numel(obj) % for each session
    % get performance for each condition
    for j = 1:numel(params(i).trialid)
        nTrialsInCond = numel(params(i).trialid{j});
        hitTrialsInCond = obj(i).bp.hit(params(i).trialid{j});
        perf(i,j) = sum(hitTrialsInCond) ./ nTrialsInCond;
    end
end

end