function data_centered = ZscoreFiringRate(data, mu, sd)
    % ZscoreFiringRate - Z-scores the firing rate data using provided mean (mu) and standard deviation (sd).
    %
    % Syntax: data_centered = ZscoreFiringRate(data, mu, sd)
    %
    % Inputs:
    %    data - The firing rate data, which can be either 2D (time x neurons) or 3D (time x trials x neurons).
    %    mu   - A vector of mean values for each neuron (neurons x 1).
    %    sd   - A vector of standard deviations for each neuron (neurons x 1).
    %
    % Outputs:
    %    data_centered - The z-scored firing rate data, same size as input data.

    % Check the dimensions of the input data
    dims = size(data);
    
    % Check if the data is 2D or 3D
    if length(dims) == 2
        % 2D data: (time x neurons)
        % Subtract the mean and divide by the standard deviation
        data_centered = (data - mu) ./ sd;
    elseif length(dims) == 3
        % 3D data: (time x trials x neurons)
        % Subtract the mean and divide by the standard deviation along the third dimension
        % Use bsxfun to apply element-wise operations along the third dimension
        data_centered = bsxfun(@minus, data, reshape(mu, [1, 1, length(mu)]));
        data_centered = bsxfun(@rdivide, data_centered, reshape(sd, [1, 1, length(sd)]));
    else
        error('Data must be either 2D (time x neurons) or 3D (time x trials x neurons)');
    end
end
