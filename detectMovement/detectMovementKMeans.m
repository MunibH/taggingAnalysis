function [moving_times, stationary_times, threshold] = detectMovementKMeans(motion_energy)
    % detectMovementKMeans - Detects when the animal is moving or stationary using K-means clustering.
    %
    % Syntax: [moving_times, stationary_times, threshold] = detectMovementKMeans(motion_energy)
    %
    % Inputs:
    %    motion_energy - A 1D array representing the motion energy time series.
    %
    % Outputs:
    %    moving_times - A binary vector of the same length as motion_energy, where 1 indicates moving and 0 indicates stationary.
    %    stationary_times - A binary vector of the same length as motion_energy, where 1 indicates stationary and 0 indicates moving.
    %    threshold - The threshold value determined using K-means clustering.
    
    % Ensure motion_energy is a column vector
    motion_energy = motion_energy(:);
    
    % Perform K-means clustering with 2 clusters (moving and stationary)
    k = 2;
    [idx, C] = kmeans(motion_energy, k);
    
    % Sort centroids to identify the stationary and moving clusters
    C_sorted = sort(C);
    stationary_centroid = C_sorted(1);
    moving_centroid = C_sorted(2);
    
    % Set the threshold as the midpoint between the two centroids
    % threshold = (stationary_centroid + moving_centroid) / 2;
    threshold = stationary_centroid;
    
    % Create binary vectors for moving and stationary periods
    moving_times = motion_energy > threshold;
    stationary_times = ~moving_times;
end
