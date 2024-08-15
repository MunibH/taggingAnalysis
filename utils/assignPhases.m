function phases = assignPhases(timeSeriesCell)
    % Input:
    % timeSeriesCell - A cell array where each cell contains a 1-D time series
    % Output:
    % phases - A cell array of the same size where each cell contains the phase
    %          values corresponding to the time points in the time series

    nSeries = length(timeSeriesCell);  % Number of time series
    phases = cell(nSeries, 1);         % Initialize output cell array
    
    for i = 1:nSeries
        nPoints = length(timeSeriesCell{i});  % Number of points in the current time series
        phases{i} = linspace(0, 2*pi, nPoints);  % Assign phases from 0 to 2*pi
    end
end
