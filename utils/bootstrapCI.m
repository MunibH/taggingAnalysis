function [ci_lower, ci_upper] = bootstrapCI(data, alpha)
    % Input:
    %   data: matrix of size (time, bootstrap iterations)
    %   alpha: confidence level (e.g., 0.05 for 95% CI)
    % Output:
    %   ci_lower: lower bound of the confidence interval (2.5th percentile)
    %   ci_upper: upper bound of the confidence interval (97.5th percentile)

    % Calculate the lower and upper percentiles
    ci_lower = prctile(data, 100 * (alpha / 2), 2);  % 2.5th percentile
    ci_upper = prctile(data, 100 * (1 - alpha / 2), 2);  % 97.5th percentile
end
