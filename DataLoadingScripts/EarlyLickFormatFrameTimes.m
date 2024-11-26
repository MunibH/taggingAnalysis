function [ft_adjusted,mask,adjustments] = EarlyLickFormatFrameTimes(ft, intervals)    
    % Remove specified time intervals from frame times and adjust the timeline
    % ft - column vector of frame times (numFrames, 1)
    % intervals - cell array of intervals, each interval is [start, end]
    % ft_adjusted - output frame times with specified intervals removed

    % Sort intervals by start time for sequential adjustments
    intervals = sortrows(cell2mat(intervals), 1);

    % Step 1: Remove frames in the specified intervals
    % Initialize mask to keep frames outside all intervals
    mask = true(size(ft));
    for i = 1:size(intervals, 1)
        mask = mask & (ft < intervals(i, 1) | ft > intervals(i, 2));
    end
    % Apply mask to remove frames
    ft_adjusted = ft(mask);

    df = diff(ft_adjusted);
    % find where df is greater than 0.05, this omits frame drops, and only
    % includes the removal points
    ix = find(df>0.005) + 1;
    adjustments = df(ix-1);

    for i = 1:numel(ix)
        ft_adjusted(ix(i):end) = ft_adjusted(ix(i):end) - adjustments(i);
    end
    % % % % % ft_adjusted(ix(2):end) = ft_adjusted(ix(2):end) - adjustments(2);
    
    % figure; plot(ft_adjusted)


end
