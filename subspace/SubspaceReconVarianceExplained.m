function [R2_per_neuron, R2_total] = SubspaceReconVarianceExplained(original_data, reconstructed_data)
    % RegressionSimilarity - Calculates the R^2 value from linear regression
    % to quantify the similarity between original and reconstructed data for each neuron and across all neurons.
    %
    % Syntax: [R2_per_neuron, R2_total] = RegressionSimilarity(original_data, reconstructed_data)
    %
    % Inputs:
    %    original_data     - A 3D array of size (time x trials x neurons) representing the original data.
    %    reconstructed_data - A 3D array of size (time x trials x neurons) representing the reconstructed data.
    %
    % Outputs:
    %    R2_per_neuron - A vector of size (neurons x 1) containing the R^2 value for each neuron.
    %    R2_total      - A scalar value representing the overall R^2 value across all neurons.

    % Check that the input dimensions match
    if ~isequal(size(original_data), size(reconstructed_data))
        error('The dimensions of the original data and reconstructed data must match.');
    end

    % Get the size of the data
    [~, ~, nNeurons] = size(original_data);

    % Preallocate the R^2 values per neuron
    R2_per_neuron = zeros(nNeurons, 1);

    % Calculate R^2 for each neuron
    for n = 1:nNeurons
        % Extract the original and reconstructed data for the current neuron
        original = original_data(:,:,n);
        reconstructed = reconstructed_data(:,:,n);

        % Reshape data into vectors
        original_vec = original(:);
        reconstructed_vec = reconstructed(:);

        % Perform linear regression: y = b0 + b1*x
        mdl = fitlm(original_vec, reconstructed_vec);

        % Extract the R^2 value from the model
        R2_per_neuron(n) = mdl.Rsquared.Ordinary;
    end

    % Calculate the overall R^2 value across all neurons
    % Combine all neurons into a single linear regression model
    original_all = original_data(:);
    reconstructed_all = reconstructed_data(:);
    mdl_total = fitlm(original_all, reconstructed_all);
    R2_total = mdl_total.Rsquared.Ordinary;
end



% function [variance_explained_per_neuron, variance_explained_total] = SubspaceReconVarianceExplained(original_data, reconstructed_data)
%     % VarianceExplained - Calculates the variance explained by the reconstructed data
%     % for each neuron and across all neurons.
%     %
%     % Syntax: [variance_explained_per_neuron, variance_explained_total] = VarianceExplained(original_data, reconstructed_data)
%     %
%     % Inputs:
%     %    original_data     - A 3D array of size (time x trials x neurons) representing the original data.
%     %    reconstructed_data - A 3D array of size (time x trials x neurons) representing the reconstructed data.
%     %
%     % Outputs:
%     %    variance_explained_per_neuron - A vector of size (neurons x 1) containing the variance explained for each neuron.
%     %    variance_explained_total      - A scalar value representing the total variance explained across all neurons.
% 
%     % Check that the input dimensions match
%     if ~isequal(size(original_data), size(reconstructed_data))
%         error('The dimensions of the original data and reconstructed data must match.');
%     end
% 
%     % Get the size of the data
%     [~, ~, nNeurons] = size(original_data);
% 
%     % Preallocate the variance explained per neuron
%     variance_explained_per_neuron = zeros(nNeurons, 1);
% 
%     % Calculate variance explained for each neuron
%     for n = 1:nNeurons
%         % Extract the original and reconstructed data for the current neuron
%         original = original_data(:,:,n);
%         reconstructed = reconstructed_data(:,:,n);
% 
%         % Calculate the variance of the original data
%         var_original = var(original(:));
% 
%         % Calculate the variance of the residual (difference between original and reconstructed)
%         residual = original - reconstructed;
%         var_residual = var(residual(:));
% 
%         % Calculate the variance explained as 1 - (variance of residual / variance of original)
%         variance_explained_per_neuron(n) = 1 - (var_residual / var_original);
%     end
% 
%     % Calculate the total variance explained across all neurons
%     total_var_original = var(original_data(:));
%     total_var_residual = var(original_data(:) - reconstructed_data(:));
%     variance_explained_total = 1 - (total_var_residual / total_var_original);
% end
