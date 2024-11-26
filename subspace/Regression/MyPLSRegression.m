function neuronWeights = MyPLSRegression(neuralData, poseTrackingData, numFolds)
    % Reshape data
    [time, neurons] = size(neuralData);
    [~, keypoints] = size(poseTrackingData);
    X = neuralData;
    Y = poseTrackingData;
    % X = reshape(neuralData, time * trials, neurons);  % Neural data
    % Y = reshape(poseTrackingData, time * trials, keypoints);  % Pose-tracking data

    % Determine maximum number of components to consider
    numComponents = min(neurons, keypoints);  

    % Initialize cross-validation
    cv = cvpartition(size(X, 1), 'KFold', numFolds);
    MSE = zeros(numComponents, keypoints);  % Track MSE for each component and keypoint

    % Manual cross-validation loop
    for k = 1:numComponents
        foldMSE = zeros(numFolds, keypoints);  % Track MSE for each fold and keypoint
        for fold = 1:numFolds
            % Split data into training and test sets
            trainIdx = cv.training(fold);
            testIdx = cv.test(fold);

            Xtrain = X(trainIdx, :);
            Ytrain = Y(trainIdx, :);
            Xtest = X(testIdx, :);
            Ytest = Y(testIdx, :);

            % Train PLS on training set with k components
            [XL, YL, XS, YS, beta] = plsregress(Xtrain, Ytrain, k);

            % Predict on test set and calculate MSE per keypoint
            Ypred = [ones(size(Xtest, 1), 1), Xtest] * beta;
            foldMSE(fold, :) = mean((Ytest - Ypred).^2, 1);  % MSE per keypoint
        end
        % Average MSE across folds for each keypoint
        MSE(k, :) = mean(foldMSE, 1);
    end

    % Select the optimal number of components for each keypoint
    [~, optimalComponentsPerKeypoint] = min(MSE);

    % Fit PLS with optimal components on entire dataset for each keypoint
    neuronWeights = zeros(neurons, keypoints);  % Initialize neuron weights
    for kp = 1:keypoints
        [XL, ~, ~, ~, ~] = plsregress(X, Y(:, kp), optimalComponentsPerKeypoint(kp));
        neuronWeights(:, kp) = sum(abs(XL(:, 1:optimalComponentsPerKeypoint(kp))), 2);  % Sum of absolute weights
    end

end