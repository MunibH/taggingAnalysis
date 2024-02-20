function dat = removeOutliers(dat,threshold)
% zscore threshold based outlier removal
% dat can be vector or matrix (time,trials) if matrix

z = abs(zscore(dat));
mask = z >= threshold;
max_dat_no_outliers = max(max(dat(~mask)));
dat(mask) = max_dat_no_outliers;

end