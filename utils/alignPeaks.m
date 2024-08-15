function alignedTimeSeries = alignPeaks(licks,activity)
    % Input:
    % timeSeries - A matrix of size (time, features)
    % Output:
    % alignedTimeSeries - A matrix of size (time, features) with peaks aligned

    [numTimePoints, numFeatures] = size(licks);
    
    % Initialize the output matrix with NaNs
    alignedTimeSeries = NaN(numTimePoints, numFeatures);

    % Find the peak position for each feature
    peakIndices = zeros(numFeatures, 1);
    for i = 1:numFeatures
        [~, peakIndices(i)] = max(licks(:, i));
    end
    
    % Determine the peak position to align to (median of all peak indices)
    peakToAlignTo = median(peakIndices);
    
    % Align each feature's peak to the chosen peak position
    for i = 1:numFeatures
        shiftAmount = peakToAlignTo - peakIndices(i);
        if shiftAmount > 0
            % Shift to the right (pad with NaNs at the beginning)
            alignedTimeSeries(1:numTimePoints-shiftAmount, i) = licks(shiftAmount+1:end, i);
        else
            % Shift to the left (pad with NaNs at the end)
            alignedTimeSeries(-shiftAmount+1:end, i) = licks(1:end+shiftAmount, i);
        end
    end
end
