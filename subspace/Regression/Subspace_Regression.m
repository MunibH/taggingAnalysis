function [in,rez] = Subspace_Regression(in,obj,input_data,input_data_zscored,right_trials,left_trials,kin)


%% Prep and Move epochs

epoch.prep(1) = mode(obj.bp.ev.goCue) - 1 - 2.15; % subtracting go cue time, 2.15
epoch.prep(2) = mode(obj.bp.ev.goCue) - 2.15;
epoch.prepix = findTimeIX(obj.time,epoch.prep,1);

epoch.move(1) = mode(obj.bp.ev.goCue) - 2.15;
epoch.move(2) = mode(obj.bp.ev.goCue) + 1 - 2.15;
epoch.moveix = findTimeIX(obj.time,epoch.move,1);

%% Get neural data

% single trial neural data
in.data.raw = input_data;
in.data.zscored = input_data_zscored;

% get trial averaged move data 
right = squeeze(mean(in.data.zscored(epoch.moveix,right_trials,:),2));
left = squeeze(mean(in.data.zscored(epoch.moveix,left_trials,:),2));
data.X = cat(1,right,left);

% perform pca on trial averaged data, keep 10 dims
[~,data.X] = pca(data.X,'NumComponents',10);

%% Get kinematic data

right = squeeze(mean(kin(epoch.moveix,right_trials,:),2));
left = squeeze(mean(kin(epoch.moveix,left_trials,:),2));
data.Y = cat(1,right,left);

%% Find W in M = WN (ridge regression)

ridge_vals = logspace(-6,6,100);
nFolds = 4;

W = MyRidgeRegression(data.X,data.Y,ridge_vals,nFolds);
% remove intercept coefficients
W = W(2:end,:);

%% Subspace ID

% primer on using svd to find spaces of a matrix
% http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf

% rank of W (how many linearly independent cols are there)
% there will be k many potent dimensions, and size(W,2)-k null dimensions
tolerance = 0.1; % rank(A,TOL) is the number of singular values of A that are larger than TOL.
k = rank(W', tolerance); 

% column, row, and null space of W can be found through SVD
[u,s,v] = svd(W'); % W' = u*s*v'. check this with the command: immse(W',u*s*v'). should return ~0

% row space of W is potent space
rez.Q.potent = v(:,1:k);

% null space of W
rez.Q.null = v(:,(k+1):end);

end











