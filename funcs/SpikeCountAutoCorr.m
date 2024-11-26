function autocorrAll = SpikeCountAutoCorr(activity, lags)
    % Calculate the spike count autocorrelation for all neurons.
    % 
    % Inputs:
    %   activity   - A (bins, trials, neurons) matrix of binned neural activity.
    %   lags       - A vector of time lags (e.g., [-10:10]).
    %
    % Output:
    %   autocorrAll - A (lags, neurons) matrix of spike count autocorrelations.
    
    % Extract dimensions
    [numBins, numTrials, numNeurons] = size(activity);
    numLags = length(lags);
    
    % Preallocate the result matrix
    autocorrAll = NaN(numLags, numNeurons);
    
    % Loop over each neuron
    for neuronIdx = 1:numNeurons
        % Extract activity for the current neuron
        neuronActivity = squeeze(activity(:, :, neuronIdx)); % (bins, trials)
        
        % Loop over each time lag
        for lagIdx = 1:numLags
            lag = lags(lagIdx);
            
            if lag >= 0
                % Positive or zero lag
                validBins = 1:(numBins - lag);
                refBins = validBins;
                laggedBins = validBins + lag;
            else
                % Negative lag
                validBins = (1 - lag):numBins;
                refBins = validBins + lag;
                laggedBins = validBins;
            end
            
            % Extract spike counts for reference and lagged bins
            refCounts = neuronActivity(refBins, :); % (validBins, trials)
            laggedCounts = neuronActivity(laggedBins, :); % (validBins, trials)
            
            % Flatten the data across bins and trials for correlation
            refCountsFlat = refCounts(:);
            laggedCountsFlat = laggedCounts(:);
            
            % Calculate correlation (skip if data is constant)
            if std(refCountsFlat) > 0 && std(laggedCountsFlat) > 0
                autocorrAll(lagIdx, neuronIdx) = corr(refCountsFlat, laggedCountsFlat);
            else
                autocorrAll(lagIdx, neuronIdx) = NaN; % Define as NaN for constant data
            end
        end
    end
end
