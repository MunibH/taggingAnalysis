function [in,rez] = Subspace_SingleTrialElsayed(in,obj,move_,sesspar,meta,params)

warning('off', 'manopt:getHessian:approx')

disp('Identifying Null and Potent Subspaces using SingleTrialElsayed for:')
disp([meta.anm ' ' meta.date ' ' meta.region])

%%

in.nProbes = numel(obj.psth);

%% Gather trials 

if strcmpi(in.trials,'all')
    trials = 1:obj.bp.Ntrials;
else
    % ensure trials is not string or char
    if isstring(in.trials) || ischar(in.trials)
        error('in.trials must be `all` or numeric conditions')
    end
    
    % get trials to use, balancing right and left trials
    trials_cond = sesspar.trialid(in.trials);
    all_trials = cell2mat(trials_cond');
    right_trials = all_trials(logical(obj.bp.R(all_trials)));
    left_trials = all_trials(logical(obj.bp.L(all_trials)));

    minTrials = min(numel(right_trials),numel(left_trials));
    
    right_trials = randsample(right_trials,minTrials,false);
    left_trials = randsample(left_trials,minTrials,false);

    trials = sort([right_trials ; left_trials]); 
end

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

%% Get neural data

if in.nProbes > 1
    % concatenate dual probes since they're same region, opposite
    % hemisphere
    input_data = cat(2,obj.trialdat{1},obj.trialdat{2}); 
    in.basemu = [obj.baseline{1}.mu ; obj.baseline{2}.mu];
    in.basesd = [obj.baseline{1}.sigma ; obj.baseline{2}.sigma];
else
    input_data = obj.trialdat{1};
    in.basemu = obj.baseline{1}.mu;
    in.basesd = obj.baseline{1}.sigma;
end
input_data = permute(input_data,[1 3 2]); % (time,trials,units)

% preprocess input_data (zscore with baseline stats)
input_data_zscored = ZscoreFiringRate(input_data,in.basemu,in.basesd);

% single trial neural data
in.data.raw = input_data;
in.data.zscored = input_data_zscored;

in.dims = size(in.data.raw);
in.data.zscored_reshaped = reshape(in.data.zscored,in.dims(1)*in.dims(2),in.dims(3));


in.data.null = in.data.zscored_reshaped(~mask,:);
in.data.potent = in.data.zscored_reshaped(mask,:);

%% Covariances

in.C.null = cov(in.data.null);
in.C.potent = cov(in.data.potent);

%% Subspace ID
% main optimization step
[Q,~,P,~,~] = orthogonal_subspaces(in.C.potent,in.nPotentDim, ...
                                   in.C.null,in.nNullDim,in.alpha);

rez.Q.potent = Q*P{1}; % potent subspace weights
rez.Q.null = Q*P{2};   % null subspace weights


end











