function [B_opt, lasso_opt, mse] = MyLassoRegression(X, Y, lasso_vals, nFolds)
    [T, DX] = size(X);
    DY = size(Y, 2);

    % Add a bias term (intercept) to X
    X = [ones(T, 1), X]; % Now X is (T, DX+1)

    % Set up cross-validation indices
    cv = cvpartition(T, 'KFold', nFolds);

    mse_lasso = zeros(length(lasso_vals), DY);
    B_opt = zeros(DX+1, DY);
    lasso_opt = zeros(1, DY);

    % Loop over each output dimension (DY)
    for dim = 1:DY
        Y_dim = Y(:, dim); % Select the current output dimension

        % Loop over each lasso parameter
        for i = 1:length(lasso_vals)
            lasso_param = lasso_vals(i);
            lasso_param = 1000;
            mse_fold = zeros(nFolds, 1);

            % Cross-validation loop
            for fold = 1:nFolds
                trainIdx = training(cv, fold);
                testIdx = test(cv, fold);

                X_train = X(trainIdx, :);
                Y_train_dim = Y_dim(trainIdx);
                X_test = X(testIdx, :);
                Y_test_dim = Y_dim(testIdx);

                % Perform lasso regression using MATLAB's lasso function
                [B, FitInfo] = lasso(X_train(:, 2:end), Y_train_dim, 'Lambda', lasso_param);

                % Include the intercept
                intercept = FitInfo.Intercept;
                B_full = [intercept; B]; 

                % Predict on test data
                Y_pred = X_test * B_full; % Now the dimensions match

                % Calculate mean squared error on test set
                mse_fold(fold) = mean((Y_test_dim - Y_pred).^2);
            end

            % Average MSE over all folds for this output dimension
            mse_lasso(i, dim) = mean(mse_fold);
        end

        % Find the optimal lasso parameter for this output dimension
        [mse_opt, opt_idx] = min(mse_lasso(:, dim));
        lasso_opt(dim) = lasso_vals(opt_idx);

        % Train on the full dataset with the optimal lasso parameter for this output dimension
        [B, FitInfo] = lasso(X(:, 2:end), Y_dim, 'Lambda', lasso_opt(dim));
        intercept = FitInfo.Intercept;
        B_opt(:, dim) = [intercept; B]; % Include intercept
    end

    mse = min(mse_lasso, [], 1);
end
