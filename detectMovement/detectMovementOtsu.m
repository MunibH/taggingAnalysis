function [moving_times, stationary_times, threshold] = detectMovementOtsu(motion_energy)
    % detectMovementOtsu - Detects when the animal is moving or stationary using Otsu's method.
    %
    % Syntax: [moving_times, stationary_times, threshold] = detectMovementOtsu(motion_energy)
    %
    % Inputs:
    %    motion_energy - A 1D array representing the motion energy time series.
    %
    % Outputs:
    %    moving_times - A binary vector of the same length as motion_energy, where 1 indicates moving and 0 indicates stationary.
    %    stationary_times - A binary vector of the same length as motion_energy, where 1 indicates stationary and 0 indicates moving.
    %    threshold - The threshold value determined using Otsu's method.
    
    % Ensure motion_energy is a column vector
    motion_energy = motion_energy(:);
    
    % Normalize motion energy to range [0, 1]
    normalized_motion_energy = (motion_energy - min(motion_energy)) / (max(motion_energy) - min(motion_energy));
    
    % Compute Otsu's threshold
    threshold = graythresh(normalized_motion_energy) * (max(motion_energy) - min(motion_energy)) + min(motion_energy);
    
    % Create binary vectors for moving and stationary periods
    moving_times = motion_energy > threshold;
    stationary_times = ~moving_times;
end