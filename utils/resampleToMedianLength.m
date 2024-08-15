function resampledData = resampleToMedianLength(timeSeriesCell)
    % Input:
    % timeSeriesCell - A cell array where each cell contains a 1-D time series
    % Output:
    % resampledData - A cell array with each time series resampled to the median length

    % Determine the number of time series
    nSeries = length(timeSeriesCell);
    
    % Get the lengths of each time series
    lengths = cellfun(@length, timeSeriesCell);
    
    % Calculate the median length
    medianLength = median(lengths);
    
    % Initialize output cell array
    resampledData = cell(nSeries, 1);
    
    % Resample each time series to the median length
    for i = 1:nSeries
        originalTime = linspace(1, medianLength, lengths(i)); % Original time points
        targetTime = 1:medianLength; % Target time points
        resampledData{i} = interp1(originalTime, timeSeriesCell{i}, targetTime, 'linear');
    end
end
