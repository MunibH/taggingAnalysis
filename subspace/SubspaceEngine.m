function SubspaceEngine(meta,inpar,obj,par,params,tag,me)

%% Gather trials

if strcmpi(inpar.trials,'all')
    all_trials = 1:obj.bp.Ntrials;
else
    % ensure trials is not string or char
    if isstring(inpar.trials) || ischar(inpar.trials)
        error('in.trials must be `all` or numeric conditions')
    end

    % get trials to use based on condition, balancing right and left trials amongst those specified trials
    trials_cond = sesspar.trialid(inpar.trials);
    all_trials = cell2mat(trials_cond');

end
% balance right and left trials
right_trials = all_trials(logical(obj.bp.R(all_trials)));
left_trials = all_trials(logical(obj.bp.L(all_trials)));

minTrials = min(numel(right_trials),numel(left_trials));

right_trials = randsample(right_trials,minTrials,false);
left_trials = randsample(left_trials,minTrials,false);

trials = sort([right_trials ; left_trials]);

%% Get neural data

[input_data,input_data_zscored] = GetSingleTrialInputData(inpar,obj);

%% SUBSPACE ID
switch inpar.method

    case 'st'
        [in,out] = Subspace_SingleTrialElsayed(inpar,obj,me.move,trials,input_data,input_data_zscored);

    case 'ta'
        [in,out] = Subspace_TrialAveragedElsayed(inpar,obj,input_data,input_data_zscored,right_trials,left_trials);

    case '2pca'
        [in,out] = Subspace_2PCA(inpar,obj,me.move,trials,input_data,input_data_zscored);

    case 'regress'
        kin = GetKinForRegression(obj,me,par);
        [in,out] = Subspace_Regression(inpar,obj,input_data,input_data_zscored,right_trials,left_trials,kin);
        
end


end











