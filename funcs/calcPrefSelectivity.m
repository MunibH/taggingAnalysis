function [sel,cluix,ix_,dprime] = calcPrefSelectivity(obj, params, cond2use, epoch, whenSelective, inp)


edges = whenSelective + mode(obj.bp.ev.(epoch));
edges = edges - mode(obj.bp.ev.(params.alignEvent));
ix = findTimeIX(obj.time,edges);
ix_ = ix(1):ix(2);

%% calculate preferred selectivity
psth = obj.psth(:,:,cond2use);

epochAvg = squeeze(mean(psth(ix_,:,:),1)); % mean over time (clu,cond)

[~,pref] = max(epochAvg'); % preferred direction
temp = pref;
nonpref = temp;
nonpref(temp==1) = 2;
nonpref(temp==2) = 1;

for i = 1:size(psth,2)
    sel(:,i) = psth(:,i,pref(i)) - psth(:,i,nonpref(i));
end

%% find which cells are significantly selective
clear epochAvg
% cluix corresponds to clu index in obj.psth/trialdat

% num trials to estimate selectivity
subTrials = inp.nTrials;
% The p-value that you want to perform the ranksum test at
sig = inp.pval;

trialdat = permute(obj.trialdat,[1 3 2]); % (time,trials,clu) 
trialdat_mean = squeeze(mean(trialdat(ix_,:,:),1)); % Take the average FR for all cells during specified time

for c = 1:numel(cond2use)
    trix = params.trialid{cond2use(c)};
    trix2use = randsample(trix,subTrials,true);
    epochAvg{c} = trialdat_mean(trix2use,:); % (trials,units)
    sampled(:,:,:,c) = trialdat(ix_,trix2use,:); % (time window,trials,units,cond)
end

nUnits = size(epochAvg{1},2);                  
pvals = zeros(1,nUnits);                     % Store p-values for each cell
hyp = zeros(1,nUnits);                       % Store hyp test results for each cell
pref = zeros(1,nUnits);
for c = 1:nUnits
    [pvals(c),hyp(c)] = ranksum(epochAvg{1}(:,c) , epochAvg{2}(:,c),'alpha',sig);
end

cluix = logical(hyp);

%% d prime


% d' num = mean_over_trials(right) - mean_over_trials(left)
% d' den = sqrt(  (var_over_trials(right) + var_over_trials(left))/2   )

n = squeeze(mean(sampled(:,:,:,1),2)) - squeeze(mean(sampled(:,:,:,2),2));
d = sqrt(  squeeze(var(sampled(:,:,:,1),[],2)) + squeeze(var(sampled(:,:,:,2),[],2))  );
dprime = n./d;








end % calcPrefSelectivity