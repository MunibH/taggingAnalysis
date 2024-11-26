function [in,rez] = Subspace_PLS_Regression(in,obj,input_data,input_data_zscored,right_trials,left_trials,kin)


%% Prep and Move epochs

epoch.prep(1) = mode(obj.bp.ev.goCue) - 1 - 2.15; % subtracting go cue time, 2.15
epoch.prep(2) = mode(obj.bp.ev.goCue) - 2.15;
epoch.prepix = findTimeIX(obj.time,epoch.prep,1);

epoch.move(1) = mode(obj.bp.ev.goCue) - 2.15;
epoch.move(2) = mode(obj.bp.ev.goCue) + 1 - 2.15;
epoch.moveix = findTimeIX(obj.time,epoch.move,1);

epoch.moveix = 1:size(input_data,1); % use all time, for debugging only

%% Get neural data

% single trial neural data
in.data.raw = input_data;
in.data.zscored = input_data_zscored;

in.C = cov(reshape(in.data.zscored,[],size(in.data.zscored,3)));

% get trial averaged move data
right = squeeze(mean(in.data.zscored(epoch.moveix,right_trials,:),2));
left = squeeze(mean(in.data.zscored(epoch.moveix,left_trials,:),2));
data.X = cat(1,right,left);

% % perform pca on trial averaged data, keep 10 dims
if in.regress.pca
    [~,data.X] = pca(data.X,'NumComponents',10);
end

%% Get kinematic data

right = squeeze(mean(kin(epoch.moveix,right_trials,:),2));
left = squeeze(mean(kin(epoch.moveix,left_trials,:),2));
data.Y = cat(1,right,left);

%% Find W using partial least squares

nFolds = 4;

W = MyPLSRegression(data.X,data.Y,nFolds);

%% Subspace ID

% primer on using svd to find spaces of a matrix
% http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf

% rank of W (how many linearly independent cols are there)
% there will be k many potent dimensions, and size(W,2)-k null dimensions
tolerance = 0.01; % rank(A,TOL) is the number of singular values of A that are larger than TOL.
k = rank(W', tolerance);

% % column, row, and null space of W can be found through SVD
% [u,s,v] = svd(W'); % W' = u*s*v'. check this with the command: immse(W',u*s*v'). should return ~0
% 
% % row space of W is potent space
% rez.Q.potent = v(:,1:k);
% 
% % null space of W
% rez.Q.null = v(:,(k+1):end);

rez.Q.potent = orth(W);
rez.Q.null = null(W');

rez.W = W;


end











