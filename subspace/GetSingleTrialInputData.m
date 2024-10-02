function [input_data,input_data_zscored] = GetSingleTrialInputData(in,obj)

switch in.standardize
    case 'baseline'
        e = 'baseline';
    case 'trial'
        e = 'trialfr';
end

if numel(obj.trialdat) > 1
    % concatenate dual probes since they're same region, opposite
    % hemisphere
    input_data = cat(2,obj.trialdat{1},obj.trialdat{2});
    in.basemu = [obj.(e){1}.mu ; obj.(e){2}.mu];
    in.basesd = [obj.(e){1}.sigma ; obj.(e){2}.sigma];
else
    input_data = obj.trialdat{1};
    in.basemu = obj.(e){1}.mu;
    in.basesd = obj.(e){1}.sigma;
end
input_data = permute(input_data,[1 3 2]); % (time,trials,units)

% preprocess input_data (zscore with baseline stats)
input_data_zscored = ZscoreFiringRate(input_data,in.basemu,in.basesd);

end