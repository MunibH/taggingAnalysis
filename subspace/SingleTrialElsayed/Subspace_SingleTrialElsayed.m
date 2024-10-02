function [in,rez] = Subspace_SingleTrialElsayed(in,obj,move_,trials,input_data,input_data_zscored)

warning('off', 'manopt:getHessian:approx')


%% Moving and stationary epochs

if in.delayOnly
    delay_t(1) = mode(obj.bp.ev.delay) - 2.5;
    delay_t(2) = mode(obj.bp.ev.goCue) - 0.02 - 2.5;
    ix = findTimeIX(obj.time,delay_t);
elseif in.responseOnly
    rt(1) = mode(obj.bp.ev.goCue) + 1 - 2.5;
    rt(2) = mode(obj.bp.ev.goCue) - 2.5;
    ix = findTimeIX(obj.time,rt);
else
    ix = 1:numel(obj.time);
end

% movement mask based off motion energy and threshold (me.move comes from loadMotionEnergy() )
mask = move_(ix,trials);
mask = mask(:); % (time*trials) , 1 where animal is moving, 0 where animal is quiet

% sample same number of move and non-move time points
nMove = sum(mask);
nNonMove = sum(~mask);
minCount = min(nMove,nNonMove);
moveix = find(mask);
nonmoveix = find(~mask);

moveix = sort(moveix(randperm(numel(moveix), minCount)));
nonmoveix = sort(nonmoveix(randperm(numel(nonmoveix), minCount)));

%% Get neural data

% single trial neural data
in.data.raw = input_data;
in.data.zscored = input_data_zscored;

in.dims = size(in.data.raw);
in.data.zscored_reshaped = reshape(in.data.zscored,in.dims(1)*in.dims(2),in.dims(3));

in.data.null = in.data.zscored_reshaped(nonmoveix,:);
in.data.potent = in.data.zscored_reshaped(moveix,:);

% calculate move and nonmove mean firing rates
input_re = reshape(input_data,in.dims(1)*in.dims(2),in.dims(3));
input_re_nonmove = input_re(nonmoveix,:);
input_re_move = input_re(moveix,:);
rez.fr.null = mean(input_re_nonmove,1)';
rez.fr.potent = mean(input_re_move,1)';

%% Covariances

in.C.null = cov(in.data.null);
in.C.potent = cov(in.data.potent);

%% Dimensionality

% % either perform parallel analysis or load results of previous parallel analysis runs 
% fn = [meta.anm '_' meta.date '_numNullPotentDims'];
% fpth = in.dimpth;
% if in.estimateDimensionality
%     disp('Estimating number of null and potent dimensions')
%     disp('Null:')
%     dims.null = ParallelAnalysis(in.data.null);
%     disp('Potent:')
%     dims.potent = ParallelAnalysis(in.data.potent);
%     if in.saveDimensionality
%         SaveResults(fpth,fn,dims);
%     end
% 
% else    
%     load(fullfile(fpth,fn))
% end
% 
% % cap dimensionality at 20
% if exist('dims','var')
%     in.nPotentDim = min(20,dims.potent);
%     in.nNullDim = min(20,dims.null);
% end


%% Subspace ID
% main optimization step
[Q,~,P,~,~] = orthogonal_subspaces(in.C.potent,in.nPotentDim, ...
                                   in.C.null,in.nNullDim,in.alpha);

rez.Q.potent = Q*P{1}; % potent subspace weights
rez.Q.null = Q*P{2};   % null subspace weights


end











