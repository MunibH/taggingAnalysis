function [sel,shadeix,selClu] = getSelectiveClu(obj, params, cond2use, epoch, whenSelective, nTrials2use)

[sel,shadeix] = calcPrefSelectivity(obj, params, cond2use, epoch, whenSelective);

edges = whenSelective + mode(obj.bp.ev.(epoch));
edges = edges - mode(obj.bp.ev.(params.alignEvent));
for i = 1:numel(whenSelective)
    [~,ix(i)] = min(abs(obj.time - edges(i)));
end
ix_ = ix(1):ix(2);

mu = mean(obj.trialdat(ix_,:,:),1);                 % (neurons,trials)Take the average FR for all cells during desired period
temp = squeeze(mu)';                                 % (trials x neurons)

for c = 1:numel(cond2use)
    trix = params.trialid{cond2use(c)};
    if numel(trix) >= nTrials2use
        rep = false;
    else
        rep = true;
    end
    trix2use = randsample(trix,nTrials2use,rep);
    
    epochAvg{c} = temp(trix,:);
end

% The p-value that you want to perform the ranksum test at
sig = 0.01;

nDims = size(epochAvg{1},2);                % Get the number of neurons
pvals = zeros(1,nDims);                     % Store p-values for each neuron
hyp = zeros(1,nDims);                       % Store hyp test results for each neuron
% pref = zeros(1,nDims);
for c = 1:nDims
    [pvals(c),hyp(c)] = ranksum(epochAvg{1}(:,c),epochAvg{2}(:,c),'alpha',sig);
%     % get preferred direction for each dimension
%     mus(1) = mean(epochAvg{1}(:,c));
%     mus(2) = mean(epochAvg{2}(:,c));
%     [~,pref(c)] = max(mus); % preferred direction is just direction with higher epochAvg
end

selClu = hyp;

end