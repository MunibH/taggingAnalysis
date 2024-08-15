function [m,h,lowerbnd,upperbnd] = mean_CI(data, confidence, varargin)
if nargin < 2
    confidence = 0.95; % Default confidence level
end

if nargin > 2
    boot = varargin{1};
else
    boot = false;
end

a = 1.0 * data;
n = size(a,2);
m = nanmean(a,2);
if boot
    den = 1;
else
    den = sqrt(n);
end
se = nanstd(a,[],2) / den;
t_critical = tinv((1 + confidence) / 2, n - 1);
h = se * t_critical;

lowerbnd = m - h;
upperbnd = m + h;



end
