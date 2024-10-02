function [in,rez] = Subspace_TrialAveragedElsayed(in,obj,input_data,input_data_zscored,right_trials,left_trials)

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

% get trial averaged preparatory data 
right = squeeze(mean(in.data.zscored(epoch.prepix,right_trials,:),2));
left = squeeze(mean(in.data.zscored(epoch.prepix,left_trials,:),2));
prepdata = cat(1,right,left);

% get trial averaged move data 
right = squeeze(mean(in.data.zscored(epoch.moveix,right_trials,:),2));
left = squeeze(mean(in.data.zscored(epoch.moveix,left_trials,:),2));
movedata = cat(1,right,left);

in.data.null = prepdata;
in.data.potent = movedata;

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











