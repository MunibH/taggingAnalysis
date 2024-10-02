function [input_data,input_data_zscored] = GetTrialAveragedInputData(in,obj)

switch in.standardize
    case 'baseline'
        e = 'baseline';
    case 'trial'
        e = 'trialfr';
end

if numel(obj.psth) > 1
    % concatenate dual probes since they're same region, opposite
    % hemisphere
    input_data = cat(2,obj.psth{1},obj.psth{2}); % (time,units,cond)
    in.basemu = [obj.(e){1}.mu ; obj.(e){2}.mu];
    in.basesd = [obj.(e){1}.sigma ; obj.(e){2}.sigma];
else
    input_data = obj.psth{1};
    in.basemu = obj.(e){1}.mu;
    in.basesd = obj.(e){1}.sigma;
end


% preprocess input_data (zscore with baseline stats)
for icond = 1:size(input_data,3)
    input_data_zscored(:,:,icond) = ZscoreFiringRate(input_data(:,:,icond),in.basemu,in.basesd);
end
input_data = permute(input_data,[1 3 2]);
input_data_zscored = permute(input_data_zscored,[1 3 2]);

end





