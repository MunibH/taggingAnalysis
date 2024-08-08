function [moving_times, stationary_times, threshold] = detectMovementGMM(motion_energy)
    % detectMovementGMM - Detects when the animal is moving or stationary using a Gaussian Mixture Model (GMM).
    %
    % Syntax: [moving_times, stationary_times, threshold] = detectMovementGMM(motion_energy)
    %
    % Inputs:
    %    motion_energy - A 1D array representing the motion energy time series.
    %
    % Outputs:
    %    moving_times - A binary vector of the same length as motion_energy, where 1 indicates moving and 0 indicates stationary.
    %    stationary_times - A binary vector of the same length as motion_energy, where 1 indicates stationary and 0 indicates moving.
    %    threshold - The threshold value determined using the GMM approach.
    
    % Ensure motion_energy is a column vector
    motion_energy = motion_energy(:);
    
    % Set GMM options
    options = statset('MaxIter', 500, 'TolFun', 1e-6, 'Display', 'final'); % Increase MaxIter and set tolerance
    
    % Fit a GMM with 2 components (stationary and moving), with regularization and specified options
    gmm = fitgmdist(motion_energy, 2, 'Options', options, 'RegularizationValue', 1e-5, 'Replicates', 5);
    
    % Extract the means of the two Gaussian components
    means = gmm.mu;
    
    % Sort the means to identify the stationary and moving components
    sorted_means = sort(means);
    stationary_mean = sorted_means(1);
    moving_mean = sorted_means(2);
    
    % Set the threshold as the midpoint between the two means
    % threshold = (stationary_mean + moving_mean) / 2;
    threshold = stationary_mean;
    
    % Create binary vectors for moving and stationary periods
    moving_times = motion_energy > threshold;
    stationary_times = ~moving_times;
end


% function [moving_times, stationary_times, threshold] = detectMovementGMM(motion_energy)
%     % detectMovementGMM - Detects when the animal is moving or stationary using a Gaussian Mixture Model (GMM).
%     %
%     % Syntax: [moving_times, stationary_times, threshold] = detectMovementGMM(motion_energy)
%     %
%     % Inputs:
%     %    motion_energy - A 1D array representing the motion energy time series.
%     %
%     % Outputs:
%     %    moving_times - A binary vector of the same length as motion_energy, where 1 indicates moving and 0 indicates stationary.
%     %    stationary_times - A binary vector of the same length as motion_energy, where 1 indicates stationary and 0 indicates moving.
%     %    threshold - The threshold value determined using the GMM approach.
% 
%     % Ensure motion_energy is a column vector
%     motion_energy = motion_energy(:);
% 
%     % Fit a GMM with 2 components (stationary and moving)
%     gmm = fitgmdist(motion_energy, 2);
% 
%     % Extract the means of the two Gaussian components
%     means = gmm.mu;
% 
%     % Sort the means to identify the stationary and moving components
%     sorted_means = sort(means);
%     stationary_mean = sorted_means(1);
%     moving_mean = sorted_means(2);
% 
%     % Set the threshold as the midpoint between the two means
%     % threshold = (stationary_mean + moving_mean) / 2;
%     threshold = stationary_mean;
% 
%     % Create binary vectors for moving and stationary periods
%     moving_times = motion_energy > threshold;
%     stationary_times = ~moving_times;
% end
