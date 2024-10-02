function SaveResults(pth, fn, varargin)
% SaveResults saves multiple variables with their original names to a .mat file.
% 
% pth  - Path to save the file
% fn   - File name
% varargin - Variables to be saved

% Check if the path exists, create it if not
if ~exist(pth, 'dir')
    mkdir(pth);
end

% Prepare the save path
savepth = fullfile(pth, [fn '.mat']);

% Loop through each input variable
for i = 1:length(varargin)
    varname = inputname(i+2); % Get the variable name (i+2 due to pth and fn)
    eval([varname, ' = varargin{i};']); % Assign the variable in the workspace
    if i==1
        save(savepth, varname, '-v7.3'); % Save the variable with its name
    else
        save(savepth, varname, '-append'); % Save the variable with its name
    end
end

disp(['Results saved to: ' savepth]);

end
