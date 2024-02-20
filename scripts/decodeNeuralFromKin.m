clear par
% DECODING PARAMETERS -----------------------------------------------------

% input data = kin data (time*trials,kinfeats)
% output data = neural data   (time*trials,units)

par.pre=15; % time bins prior to output used for decoding
par.post=2; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

% data sets
par.train = 0.7; % fraction of trials
par.test = 1 - par.train;

% feature to use to decode
par.feats = kin(1).featLeg;
par.feats = {'motion','nose','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});

% trials
par.cond2use = [2, 3];
par.conds = {'rhit','lhit'};

% DECODING -----------------------------------------------------

% trials
for i = 1:numel(par.cond2use)
    par.trials.(par.conds{i}) = params.trialid{par.cond2use(i)};
    nTrialsCond(i) = numel(par.trials.(par.conds{i}));
end
% balance number of trials across conditions
minTrials = min(nTrialsCond);
mask = nTrialsCond > minTrials;
par.trials.all = [];
for i = 1:numel(par.cond2use)
    if mask(i)
        par.trials.(par.conds{i}) = sort(randsample(par.trials.(par.conds{i}), minTrials, false),'ascend');
    end
    par.trials.all = [par.trials.all; par.trials.(par.conds{i})];
end

% partition train and test
nTrials = numel(par.trials.all);
nTrain = floor(nTrials*par.train);
par.trials.train = randsample(par.trials.all,nTrain,false);
par.trials.test = par.trials.all(~ismember(par.trials.all,par.trials.train));

% input data
par.featix = find(ismember(kin.featLeg,par.feats));

X.train = kin.dat(:,par.trials.train,par.featix); % (time,trials,feats)
X.size.train = size(X.train);
X.train = reshape(X.train, size(X.train,1)*size(X.train,2),size(X.train,3));

X.test = kin.dat(:,par.trials.test,par.featix); % (time,trials,feats)
X.size.test = size(X.test);
X.test = reshape(X.test, size(X.test,1)*size(X.test,2),size(X.test,3));

% reshape train and test data to account for prediction bin size
X.train = reshapePredictors(X.train,par);
X.test = reshapePredictors(X.test,par);

% flatten inputs
% if you're using a model with recurrence, don't flatten
X.train = reshape(X.train,size(X.train,1),size(X.train,2)*size(X.train,3));
X.test = reshape(X.test,size(X.test,1),size(X.test,2)*size(X.test,3));

% output data
Y.train = permute(obj.trialdat(:,:,par.trials.train), [1 3 2]); % (time,trials,units);
Y.size.train = size(Y.train);
Y.train = reshape(Y.train, size(Y.train,1)*size(Y.train,2),size(Y.train,3));

Y.test = permute(obj.trialdat(:,:,par.trials.test), [1 3 2]); % (time,trials,units);
Y.size.test = size(Y.test);
Y.test = reshape(Y.test, size(Y.test,1)*size(Y.test,2),size(Y.test,3));

% standardize data
% standardize both train and test sets using train set statistics
% can also standardize using specific time points (presample for example)
X.mu = mean(X.train,1,'omitnan');
X.sigma = std(X.train,[],1,'omitnan');
X.train = (X.train - X.mu) ./ X.sigma;
if ~par.test==0
    X.test = (X.test - X.mu) ./ X.sigma;
end

Y.mu = mean(Y.train,1,'omitnan');
Y.sigma = std(Y.train,[],1,'omitnan');
Y.train = (Y.train - Y.mu) ./ Y.sigma;
if ~par.test==0
    Y.test = (Y.test - Y.mu) ./ Y.sigma;
end

% fill missing values in kinematics
X.train = fillmissing(X.train,'constant',0);
Y.train = fillmissing(Y.train,'nearest');
X.test = fillmissing(X.test,'constant',0);
Y.test = fillmissing(Y.test,'nearest');

% train and test models
for i = 1:Y.size.train(3)
    if mod(i,20) == 0
        disp(['Unit ' num2str(i) '/' num2str(Y.size.train(3))])
    end
    mdl{i} = fitrlinear(X.train, Y.train(:,i), 'Learner', 'leastsquares','Regularization','ridge');
    % mdl{i} = fitrlinear(X.train,Y.train(:,i));
    Y.pred(:,i) = predict(mdl{i},X.test);
end

% variance explained
y = reshape(Y.test,Y.size.test(1),Y.size.test(2),[]); % original input data (centered)
yhat = reshape(Y.pred,Y.size.test(1),Y.size.test(2),[]); % prediction
for i = 1:size(y,3)
    tempcorr = corrcoef(y(:,i),yhat(:,i));
    Y.R2(i) = tempcorr(1,2).^2;
end

%% predict activity for all trials
clear XX YY

XX = kin.dat(:,:,par.featix); % (time,trials,feats)
XX = reshape(XX, size(XX,1)*size(XX,2),size(XX,3));

% reshape train and test data to account for prediction bin size
XX = reshapePredictors(XX,par);

% flatten inputs
% if you're using a model with recurrence, don't flatten
XX = reshape(XX,size(XX,1),size(XX,2)*size(XX,3));

% output data
YY.data = permute(obj.trialdat, [1 3 2]); % (time,trials,units);
YY.data = reshape(YY.data, size(YY.data,1)*size(YY.data,2),size(YY.data,3));

% standardize data
% standardize both train and test sets using train set statistics
% can also standardize using specific time points (presample for example)
XXmu = mean(XX,1,'omitnan');
XXsigma = std(XX,[],1,'omitnan');
XX = (XX - XXmu) ./ XXsigma;

YY.mu = mean(YY.data,1,'omitnan');
YY.sigma = std(YY.data,[],1,'omitnan');
YY.data = (YY.data - YY.mu) ./ YY.sigma;

% fill missing values in kinematics
XX = fillmissing(XX,'constant',0);

% predictions
for i = 1:numel(mdl)
    if mod(i,20) == 0
        disp(['Unit ' num2str(i) '/' num2str(Y.size.train(3))])
    end
    YY.pred(:,i) = predict(mdl{i},XX);
end

y = reshape(YY.data,numel(obj.time),obj.bp.Ntrials,[]); % original input data (centered)
yhat = reshape(YY.pred,numel(obj.time),obj.bp.Ntrials,[]); % prediction (time,trials,units)

% variance explained
for i = 1:size(y,3)
    tempcorr = corrcoef(y(:,i),yhat(:,i));
    YY.R2(i) = tempcorr(1,2).^2;
end