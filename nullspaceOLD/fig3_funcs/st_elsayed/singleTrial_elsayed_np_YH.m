function rez = singleTrial_elsayed_np_YH(input_data,obj,me,params,cond2use, varargin)

warning('off', 'manopt:getHessian:approx')

if nargin > 5
    nNPDims = varargin{1};
else
    nNPDims = nan;
end

%% trials to use

trials_cond = params.trialid(cond2use);
trials = trials_cond{1}; % all trials

%% split data into quiet and moving time points


% motion energy
mask = me.move(:,trials);
mask = mask(:); % (time*trials) , 1 where animal is moving, 0 where animal is quiet

% single trial neural data
N.full = input_data(:,trials,:);
N.dims = size(N.full);
N.full_reshape = reshape(N.full,N.dims(1)*N.dims(2),N.dims(3));


N.null = N.full_reshape(~mask,:);

mask_ = mask;
N.potent = N.full_reshape(mask_,:);

rez.N = N;

%% null and potent spaces

% -----------------------------------------------------------------------
% -- compute covariance matrices --
% -----------------------------------------------------------------------

% % method 1 - recover covariance estimate from factor analysis
% [lambda,psi] = factoran(N.null,10);
% rez.covNull = lambda*lambda' + psi;
% [lambda,psi] = factoran(N.potent,10);
% rez.covPotent = lambda*lambda' + psi;

% % method 2 - standard method
rez.covNull = cov(N.null);
rez.covPotent = cov(N.potent);


% -----------------------------------------------------------------------
% -- number of null and potent dims --
% -----------------------------------------------------------------------
% assign num dims by amount of PCs needed to explain some amount of
% variance defined in params (capped at 10 dims)
rez.varToExplain = 80;
% [rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
% if rez.dPrep > 10; rez.dPrep = 10; end
% if rez.dMove > 10; rez.dMove = 10; end
% if rez.dPrep == 1; rez.dPrep = 10; end
% if rez.dMove == 1; rez.dMove = 10; end
if isnan(nNPDims)
    rez.dPrep = floor(size(rez.covNull,1)/2);
    rez.dMove = ceil(size(rez.covNull,1)/2);
    if rez.dPrep > 20; rez.dPrep = 20; end
    if rez.dMove > 20; rez.dMove = 20; end
else
    rez.dPrep = nNPDims;
    rez.dMove = nNPDims;
end
% method 2 - keep reducing var2explain until dMove+dPrep <= full dim
% check = 1;
% while check
%     [rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
%     if (rez.dPrep + rez.dMove) <= size(N.full_reshape,2)
%         break
%     end
%     rez.varToExplain = rez.varToExplain - 1;
% end




% main optimization step
alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.covPotent,rez.dMove, ...
    rez.covNull,rez.dPrep,alpha);

rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

end











