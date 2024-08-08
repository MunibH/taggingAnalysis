function yout = decodeNeuralFromKin_v2(obj,meta,params,par,trials,alltrials,kindat,tagdat)
% DECODING -----------------------------------------------------

% only using a train set right now, because I want to predict and subtract
% predictions for all trials

% input data
X.train = kindat;
X.size.train = size(X.train);
X.train = reshape(X.train,X.size.train(1)*X.size.train(2),X.size.train(3));

% reshape train and test data to account for prediction bin size
X.train = reshapePredictors(X.train,par);

% flatten inputs
% if you're using a model with recurrence, don't flatten
X.train = reshape(X.train,size(X.train,1),size(X.train,2)*size(X.train,3));

% output data
Y.train = tagdat;
Y.size.train = size(Y.train);
Y.train = reshape(Y.train, size(Y.train,1)*size(Y.train,2),size(Y.train,3));

% standardize data
% standardize both train and test sets using train set statistics
% can also standardize using specific time points (presample for example)
X.mu = mean(X.train,1,'omitnan');
X.sigma = std(X.train,[],1,'omitnan');
X.train = (X.train - X.mu) ./ X.sigma;

Y.mu = mean(Y.train,1,'omitnan');
Y.sigma = std(Y.train,[],1,'omitnan');
Y.train = (Y.train - Y.mu) ./ Y.sigma;

% fill missing values in kinematics
X.train = fillmissing(X.train,'constant',0);
Y.train = fillmissing(Y.train,'nearest');

% train and test models
mdl = fitrlinear(X.train, Y.train, 'Learner', 'leastsquares','Regularization','ridge');
Y.pred = predict(mdl,X.train);

Y.pred = (Y.pred + Y.mu) .* Y.sigma;
Y.train = (Y.train + Y.mu) .* Y.sigma;

% variance explained
y = reshape(Y.train,Y.size.train(1),Y.size.train(2),[]); % original input data 
yhat = reshape(Y.pred,Y.size.train(1),Y.size.train(2),[]); % prediction
tempcorr = corrcoef(y(:),yhat(:));
Y.R2 = tempcorr(1,2).^2;

yout.input = y;
yout.pred = yhat;
yout.r2 = Y.R2;

end