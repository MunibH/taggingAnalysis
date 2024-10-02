function [B_opt, ridge_opt, mse] = MyRidgeRegression(X, Y, ridge_vals, nFolds)
    % Inputs:
    % X - (T, DX) input time series
    % Y - (T, DY) output time series
    % ridge_vals - vector of ridge parameters to test (e.g., logspace(-6, 6, 100))
    % nFolds - number of cross-validation folds

    % Outputs:
    % B_opt - optimal regression coefficients (DX+1, DY) including intercept
    % ridge_opt - optimal ridge regularization parameter
    % mse - mean squared error of the prediction using the optimal ridge parameter

    [T, DX] = size(X);
    DY = size(Y, 2);

    % Add a bias term (intercept) to X by concatenating a column of ones
    X = [ones(T, 1), X]; % Now X is (T, DX+1)

    % Set up cross-validation indices
    cv = cvpartition(T, 'KFold', nFolds);

    mse_ridge = zeros(length(ridge_vals), 1); % Mean squared errors for each ridge value

    % Loop over each ridge parameter
    for i = 1:length(ridge_vals)
        ridge_param = ridge_vals(i);
        mse_fold = zeros(nFolds, 1); % Store MSE for each fold

        % Cross-validation loop
        for fold = 1:nFolds
            trainIdx = training(cv, fold);
            testIdx = test(cv, fold);

            X_train = X(trainIdx, :);
            Y_train = Y(trainIdx, :);
            X_test = X(testIdx, :);
            Y_test = Y(testIdx, :);

            % Perform ridge regression: B = inv(X'X + lambda*I) X'Y
            B = (X_train' * X_train + ridge_param * eye(DX+1)) \ (X_train' * Y_train);

            % Predict on test data
            Y_pred = X_test * B;

            % Calculate mean squared error on test set
            mse_fold(fold) = mean(mean((Y_test - Y_pred).^2));
        end

        % Average MSE over all folds
        mse_ridge(i) = mean(mse_fold);
    end

    % Find the optimal ridge parameter
    [mse, opt_idx] = min(mse_ridge);
    ridge_opt = ridge_vals(opt_idx);

    % Train on the full dataset with the optimal ridge parameter
    B_opt = (X' * X + ridge_opt * eye(DX+1)) \ (X' * Y);
end
