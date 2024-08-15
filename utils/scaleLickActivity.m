function scaledActivity = scaleLickActivity(lickts, lick_duration, activity, targetDuration)
    % Input:
    % lickts - Time series of size (time, nLicks)
    % lick_duration - Vector of size (nLicks, 1) containing durations of each lick
    % activity - Time series of size (time, nLicks)
    % targetDuration - The duration to which each lick's activity will be scaled
    % Output:
    % scaledActivity - Scaled time series of size (targetDuration, nLicks)

    [numTimePoints, nLicks] = size(activity);
    
    % Initialize the output matrix
    scaledActivity = NaN(targetDuration, nLicks);
    
    % Time vectors for interpolation
    originalTime = (1:numTimePoints)';
    targetTime = linspace(1, numTimePoints, targetDuration)';
    
    for i = 1:nLicks
        % Get the duration of the current lick
        duration = lick_duration(i);
        
        % Interpolate the activity data for the current lick
        if duration < targetDuration
            % If duration is less than targetDuration, pad with NaNs at the end
            scaledActivity(1:duration, i) = interp1(originalTime, activity(:, i), linspace(1, numTimePoints, duration), 'linear', 'extrap');
        else
            % If duration is more than or equal to targetDuration, interpolate to the targetDuration
            scaledActivity(:, i) = interp1(originalTime, activity(:, i), targetTime, 'linear', 'extrap');
        end
    end
end