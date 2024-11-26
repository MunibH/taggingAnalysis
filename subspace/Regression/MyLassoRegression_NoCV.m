function [B_opt, mse] = MyLassoRegression_NoCV(X, Y)
    [T, DX] = size(X);
    DY = size(Y, 2);
    lasso_param = 1e-2; % Hardcoded lasso parameter

    % Add a bias term (intercept) to X
    X = [ones(T, 1), X]; % Now X is (T, DX+1)

    % Initialize output variables
    B_opt = zeros(DX+1, DY);
    mse = zeros(1, DY);

    % Loop over each output dimension (DY)
    for dim = 1:DY
        Y_dim = Y(:, dim); % Select the current output dimension

        % Perform lasso regression using MATLAB's lasso function
        [B, FitInfo] = lasso(X(:, 2:end), Y_dim, 'Lambda', lasso_param);

        % Include the intercept
        intercept = FitInfo.Intercept;
        B_opt(:, dim) = [intercept; B]; % Include intercept

        % Calculate mean squared error for the whole dataset
        Y_pred = X * B_opt(:, dim); % Predicted Y
        mse(dim) = mean((Y_dim - Y_pred).^2); % Mean squared error
    end
end
