function ve = SubspaceReconVarianceExplained(data, recon)
% This function calculates the variance explained for each neuron
% Inputs:
%   data     - actual data matrix of size (time, trials, neurons)
%   data_hat - model predicted data matrix of the same size (time, trials, neurons)
% Output:
%   var_explained - variance explained for each neuron (1 x neurons)

% Calculate the total variance for each neuron
total_var = squeeze(var(reshape(data, [], size(data, 3)), 0, 1)); % Variance across time and trials

if isstruct(recon)
    fns = fieldnames(recon);
    for i = 1:numel(fns)
    
        % Calculate the error variance for each neuron (difference between actual and predicted)
        residuals = data - recon.(fns{i});
        error_var = squeeze(var(reshape(residuals, [], size(data, 3)), 0, 1));
    
        % Variance explained for each neuron
        vetemp = 1 - (error_var ./ total_var);
    
        % Handle edge cases where total variance is zero
        vetemp(total_var == 0) = NaN;
        ve.(fns{i}) = vetemp;
    end
   
else
    % Calculate the error variance for each neuron (difference between actual and predicted)
    residuals = data - recon;
    error_var = squeeze(var(reshape(residuals, [], size(data, 3)), 0, 1));

    % Variance explained for each neuron
    ve = 1 - (error_var ./ total_var);

    % Handle edge cases where total variance is zero
    ve(total_var == 0) = NaN;

end
end

