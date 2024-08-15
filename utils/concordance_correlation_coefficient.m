function ccc = concordance_correlation_coefficient(x, y)
    % Calculate means
    mean_x = mean(x);
    mean_y = mean(y);
    
    % Calculate variances
    var_x = var(x);
    var_y = var(y);
    
    % Calculate covariance
    cov_xy = cov(x, y);
    cov_xy = cov_xy(1, 2);
    
    % Calculate Pearson correlation
    rho = cov_xy / (std(x) * std(y));
    
    % Calculate CCC
    ccc = (2 * rho * std(x) * std(y)) / (var_x + var_y + (mean_x - mean_y)^2);
end
