function [cd,data,par] = calcCDAll(obj,params,par)
%%

% trials
for i = 1:numel(par.cond2use)
    par.trials.cond{i} = params.trialid{par.cond2use(i)};
    nTrialsCond(i) = numel(par.trials.cond{i});
end
% balance number of trials across conditions
minTrials = min(nTrialsCond);
nTrain = floor(minTrials * 0.7);
nTest = minTrials - nTrain;

% mask = nTrialsCond > minTrials;
% par.trials.all = [];
% for i = 1:numel(par.cond2use)
%     if mask(i)
%         par.trials.cond{i} = sort(randsample(par.trials.cond{i}, minTrials, false),'ascend');
%     end
%     par.trials.all = [par.trials.all; par.trials.cond{i}];
% end

% partition train and test
for i = 1:numel(par.cond2use)
    par.trials.train{i} = randsample(par.trials.cond{i},nTrain,false);
    par.trials.test{i} = randsample(par.trials.cond{i}(~ismember(par.trials.cond{i},par.trials.train{i})),nTest,false);
end
par.trials.alltrain = cell2mat(par.trials.train');


% get train and test data
data.train = cellfun(@(x) obj.trialdat(:,:,x),par.trials.train,'uni',0);
data.test = cellfun(@(x) obj.trialdat(:,:,x),par.trials.test,'uni',0);

% zscore data
mu = nanmean( cat(3,data.train{:}),3 );
std_ = nanstd( cat(3,data.train{:}),[],3 );
data.ztrain = cellfun(@(x) (x-mu), data.train, 'uni', 0);
% data.ztrain = cellfun(@(x) (x-mu)./std_, data.train, 'uni', 0);

mu = nanmean( cat(3,data.test{:}),3 );
std_ = nanstd( cat(3,data.test{:}),[],3 );
data.ztest = cellfun(@(x) (x-mu), data.test, 'uni', 0);
% data.ztest = cellfun(@(x) (x-mu)./std_, data.test, 'uni', 0);

% calculate PSTHs
for i = 1:numel(par.cond2use)
    data.psth.train(:,:,i) = nanmean(data.ztrain{i},3);
    data.psth.test(:,:,i) = nanmean(data.ztest{i},3);
end

% calculate CD at each time point
a = data.psth.train(:,:,1) - data.psth.train(:,:,2); % (time,units)
for i = 1:numel(obj.time)
    n(i) = norm(a(i,:));
end
n = repmat(n',1,size(a,2));

cd = a./n;
cd = cd ./ norm(cd,1);

%% projections

for itime = 1:numel(obj.time)
    for i = 1:numel(par.cond2use)
        data.proj(itime,:,i) = squeeze(data.ztest{i}(itime,:,:))' * cd(itime,:)'; % (time,trials)
    end
end


%%
% 
% f = figure;
% ax = gca;
% hold on;
% plot(obj.time,data.psth.test(:,:,1),'Color','b')
% plot(obj.time,data.psth.test(:,:,2),'Color','r')

end % calcCDAll





















