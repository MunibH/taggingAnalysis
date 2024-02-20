function renummed_startAt1 = renum(input_array)
%     # https://stackoverflow.com/questions/62169627/renumber-a-sequence-to-remove-gaps-but-keep-identical-numbers
%     # takes an array that has gaps in a sequence and removes gaps, but keeps duplicates.  
% #     renum(np.array([1, 1, 1, 2, 2, 2]))  # already in correct shape
% #     > [1, 1, 1, 2, 2, 2]
% 
% #     renum(np.array([1, 1, 2, 2, 4, 4, 5, 5, 5]))  # A jump between 2 and 4
% #     > [1,1, 2, 2, 3, 3, 4, 4, 4]
% 
% #     renum(np.array([1, 1, 2, 2, 5, 2, 2]))  # A forward and backward jump
% #     > [1,1, 2, 2, 3, 4, 4]

    % make column vector if not,
    if size(input_array,1)>size(input_array,2)
        input_array = input_array';
    end
    diff_vals = diff(double(input_array));
    diff_vals(diff_vals ~= 0) = 1;

    csum_diff_vals = cumsum(diff_vals);
    if csum_diff_vals(1)<input_array(1)
        csum_diff_vals = csum_diff_vals + input_array(1);
    end
    renummed = [input_array(1), csum_diff_vals]; 
    renummed_startAt1 = renummed - (renummed(1) - 1);
end