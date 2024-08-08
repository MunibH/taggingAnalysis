function zeroCrossings = findZeroCrossings(timeSeries)
    % findZeroCrossings identifies the indices of zero crossings in a 1-D time series
    %
    % Input:
    %   timeSeries - a 1-D array of numerical values
    %
    % Output:
    %   zeroCrossings - indices where zero crossings occur

    % Ensure the input is a column vector
    if isrow(timeSeries)
        timeSeries = timeSeries';
    end
    
    % Compute the sign of each element in the time series
    signSeries = sign(timeSeries);
    
    % Find the indices where the sign changes
    zeroCrossingIndices = find(diff(signSeries) ~= 0);
    
    % Output the zero crossing indices
    zeroCrossings = zeroCrossingIndices;
end