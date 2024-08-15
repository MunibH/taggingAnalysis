function var_explained_binned = variance_explained_binned(y_true, y_pred, bin_size)
    % Calculate the variance explained in overlapping bins over time.
    %
    % Parameters:
    % - y_true: matrix of size (time, features) with original data.
    % - y_pred: matrix of size (time, features) with predicted data.
    % - bin_size: integer, the number of time points per bin.
    % - overlap: integer, the number of overlapping time points between consecutive bins.
    %
    % Returns:
    % - var_explained_binned: matrix of size (n_bins, features) with variance explained for each bin.

    [time_points,~, num_units] = size(y_true);    
    

    for i = time_points:-1:1
        end_idx = i;
        start_idx = max(end_idx - bin_size + 1, 1); % Ensure end_idx does not exceed bounds

        y_true_bin = y_true(start_idx:end_idx, :, :);
        y_true_bin = reshape(y_true_bin,size(y_true_bin,1)*size(y_true_bin,2),size(y_true_bin,3)); % (bins*trials,features)
        y_pred_bin = y_pred(start_idx:end_idx, :, :);
        y_pred_bin = reshape(y_pred_bin,size(y_pred_bin,1)*size(y_pred_bin,2),size(y_pred_bin,3)); % (bins*trials,features)

        ss_res = mean((y_true_bin - y_pred_bin) .^ 2, 1);
        ss_tot = mean((y_true_bin - mean(y_true_bin, 1)) .^ 2, 1);

        % var_explained_binned(i, :) = 1 - ss_res ./ ss_tot;
        
        for j = 1:num_units
            temp = corrcoef(y_true_bin(:,j),y_pred_bin(:,j));
            var_explained_binned(i,j) = temp(1,2).^2;
        end
    end
end

% % function var_explained_binned = variance_explained_binned(y_true, y_pred, bin_size)
% %     % Calculate the variance explained in bins over time.
% %     %
% %     % Parameters:
% %     % - y_true: matrix of size (time,trials, features) with original data.
% %     % - y_pred: matrix of size (time,trials, features) with predicted data.
% %     % - bin_size: integer, the number of time points per bin.
% %     %
% %     % Returns:
% %     % - var_explained_binned: matrix of size (n_bins, features) with variance explained for each bin.
% % 
% %     [time_points,~, num_units] = size(y_true);
% %     n_bins = floor(time_points / bin_size);
% % 
% %     var_explained_binned = zeros(n_bins, num_units);
% % 
% %     for i = 1:n_bins
% %         start_idx = (i-1) * bin_size + 1;
% %         end_idx = start_idx + bin_size - 1;
% % 
% %         y_true_bin = y_true(start_idx:end_idx, :, :);
% %         y_true_bin = reshape(y_true_bin,size(y_true_bin,1)*size(y_true_bin,2),size(y_true_bin,3)); % (bins*trials,features)
% %         y_pred_bin = y_pred(start_idx:end_idx, :, :);
% %         y_pred_bin = reshape(y_pred_bin,size(y_pred_bin,1)*size(y_pred_bin,2),size(y_pred_bin,3)); % (bins*trials,features)
% % 
% %         % ss_res = sum((y_true_bin - y_pred_bin) .^ 2, 1);
% %         % ss_tot = sum((y_true_bin - mean(y_true_bin, 1)) .^ 2, 1);
% %         % 
% %         % var_explained_binned(i, :) = 1 - ss_res ./ ss_tot;
% % 
% %         for j = 1:num_units
% %             temp = corrcoef(y_true_bin(:,j),y_pred_bin(:,j));
% %             var_explained_binned(i,j) = temp(1,2).^2;
% %         end
% %     end
% % end
