function med = findMedianJawTimes(med,jawStart,nLicks)

med.lickStart = nan(nLicks,1);

for i = 1:nLicks
    % find trials with atleast (i-1) licks
    trix_w_licks = cell2mat(cellfun(@(x) (numel(x)>(i-1)), jawStart, 'UniformOutput', false)); 
    % get i_th lick start, end, and duration for each of trix_w_licks
    ith_lickStart = cell2mat(cellfun(@(x) x(i), jawStart(trix_w_licks), 'UniformOutput', false));
    % calculate median time for i_th lick
    med.lickStart(i) = median(ith_lickStart);
end

end
