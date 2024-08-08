function restingPosition = findRestingJawPosition(timeSeries)
    % findRestingJawPosition identifies the resting position of the jaw
    % where the signal is not changing for the longest amount of time
    %
    % Input:
    %   timeSeries - a 1-D array of numerical values representing the jaw position
    %
    % Output:
    %   restingPosition - the value of the jaw position where it remains constant the longest

    % Ensure the input is a column vector
    if isrow(timeSeries)
        timeSeries = timeSeries';
    end
    
    % Tolerance for detecting flat segments (change threshold)
    tolerance = 1e-3;
    
    % Find the differences between consecutive elements
    diffSeries = abs(diff(timeSeries));
    
    % Identify where the differences are below the tolerance
    flatSegments = diffSeries < tolerance;
    
    % Append false at the beginning and end for segment detection
    flatSegments = [false; flatSegments; false];
    
    % Identify the start and end of flat segments
    segmentStarts = find(diff(flatSegments) == 1);
    segmentEnds = find(diff(flatSegments) == -1) - 1;
    
    % Calculate the lengths of the flat segments
    segmentLengths = segmentEnds - segmentStarts + 1;
    
    % Find the longest flat segment
    [~, longestSegmentIndex] = max(segmentLengths);
    
    % Determine the resting position as the mean value of the longest flat segment
    if isempty(longestSegmentIndex)
        restingPosition = NaN;
        warning('No flat segment found in the time series.');
    else
        longestSegmentStart = segmentStarts(longestSegmentIndex);
        longestSegmentEnd = segmentEnds(longestSegmentIndex);
        restingPosition = mean(timeSeries(longestSegmentStart:longestSegmentEnd));
    end
end
