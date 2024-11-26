function [accuracy, predictions] = FacemapDecodeTrialType(trialTypes, neuralActivity, windowSize)
% Decode trial type from neural activity separately at each time point
% using logistic regression with 4-fold cross-validation.
%
% INPUTS:
% trialTypes      - Vector of trial types coded as -1 or 1, size (trials, 1)
% neuralActivity  - 3D matrix of neural activity, size (time, trials, neurons)
%
% OUTPUTS:
% accuracy        - Vector of classification accuracies for each time point (size: time, 1)
% predictions     - Matrix of predicted trial types for each trial at each time point (size: time, trials)

% Get dimensions
[timePoints, numTrials, numNeurons] = size(neuralActivity);

% Ensure trialTypes is a column vector
trialTypes = trialTypes(:);

% Initialize k-fold cross-validation (k=4)
k = 4;
cv = cvpartition(numTrials, 'KFold', k);

% Preallocate results
accuracy = zeros(timePoints, 1);         % Store accuracy for each time point
predictions = NaN(timePoints, numTrials); % Store predictions for each trial and time point

% Loop through each time point
for t = 1:timePoints
    % Define the time window around the current time point t
    startTime = max(1, t - windowSize);   % Ensure window doesn't go below index 1
    endTime = min(timePoints, t + windowSize); % Ensure window doesn't exceed the total time points

    % Extract neural activity for this time point (trials x neurons)
    X = squeeze(neuralActivity(startTime:endTime, :, :))'; % Size: (trials, neurons)

    if any(all(isnan(X)))
        continue
    end

    % Store accuracy for each fold
    foldAccuracies = zeros(k, 1);

    % Store predictions for each trial in this time point
    timePointPredictions = zeros(numTrials, 1);

    % Loop over each fold for cross-validation
    for fold = 1:k
        % Get training and testing indices
        trainIdx = training(cv, fold);
        testIdx = test(cv, fold);

        % Training data and labels
        Xtrain = X(trainIdx, :);
        ytrain = trialTypes(trainIdx);

        % Testing data and labels
        Xtest = X(testIdx, :);
        ytest = trialTypes(testIdx);

        % Train logistic regression classifier
        mdl = fitclinear(Xtrain, ytrain, 'Learner', 'logistic');

        % Predict the trial types for the test data
        ypred = predict(mdl, Xtest);

        % Store predictions
        timePointPredictions(testIdx) = ypred;

        % Calculate accuracy for this fold
        foldAccuracies(fold) = mean(ypred == ytest);
    end

    % Compute overall accuracy for this time point (mean of fold accuracies)
    accuracy(t) = mean(foldAccuracies);

    % Store the predictions for this time point
    predictions(t, :) = timePointPredictions;

end
end
