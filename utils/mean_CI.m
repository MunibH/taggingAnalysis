function [m,h,lowerbnd,upperbnd] = mean_CI(data, confidence)
if nargin < 2
    confidence = 0.95; % Default confidence level
end

a = 1.0 * data;
n = size(a,2);
m = mean(a,2);
se = std(a,[],2) / sqrt(n);
t_critical = tinv((1 + confidence) / 2, n - 1);
h = se * t_critical;

lowerbnd = m - h;
upperbnd = m + h;



end
