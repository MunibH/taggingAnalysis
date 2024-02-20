function dprime = calcdprime(obj, params, cond2use, inp)

% output - dprime = (time,units)

%% d prime

% num trials to estimate selectivity
subTrials = inp.nTrials;

trialdat = permute(obj.trialdat,[1 3 2]); % (time,trials,clu) 

for c = 1:numel(cond2use)
    trix = params.trialid{cond2use(c)};
    trix2use = randsample(trix,subTrials,true);
    sampled(:,:,:,c) = trialdat(:,trix2use,:); % (time,trials,units,cond)
end

% d' num = mean_over_trials(right) - mean_over_trials(left)
% d' den = sqrt(  (var_over_trials(right) + var_over_trials(left))/2   )

n = squeeze(mean(sampled(:,:,:,1),2)) - squeeze(mean(sampled(:,:,:,2),2));
d = sqrt(  squeeze(var(sampled(:,:,:,1),[],2)) + squeeze(var(sampled(:,:,:,2),[],2))  );
dprime = (n./d);

end % calcdprime