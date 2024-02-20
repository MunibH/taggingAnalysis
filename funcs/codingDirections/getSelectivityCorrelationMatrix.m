function sel_corr_mat = getSelectivityCorrelationMatrix(sel)

% sel (time, clu)

sel_corr_mat = zeros(size(sel,1),size(sel,1));

for i = 1:size(sel_corr_mat,1) % time
    for j = 1:size(sel_corr_mat,1) % time
        temp = corrcoef(sel(i,:),sel(j,:));
        sel_corr_mat(i,j) = temp(1,2);
    end
end


end