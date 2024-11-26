
function mostRecentFile = findMostRecentFile(directory, tokens)
% Ensure tokens is a cell array of strings
if ischar(tokens) || isstring(tokens)
    tokens = {tokens}; % Convert single token to cell array
elseif ~iscell(tokens)
    error('Tokens must be a cell array, string, or character array.');
end

% Get a list of all files in the directory
files = dir(directory);

% Filter files based on the tokens
matchingFiles = files(arrayfun(@(file) all(cellfun(@(token) contains(file.name, token), tokens)), files));

if isempty(matchingFiles)
    error('No files found containing all specified tokens.');
end

% Extract the dates from the filenames
fileDates = zeros(length(matchingFiles), 1);
for i = 1:length(matchingFiles)
    % Assume the date is at the end of the filename before the extension
    [~, name, ~] = fileparts(matchingFiles(i).name);

    % Find the last set of 8 consecutive digits (assumed to be yyyymmdd)
    dateToken = regexp(name, '\d{8}$', 'match');

    if ~isempty(dateToken)
        % Convert the dateToken string to a number for comparison
        fileDates(i) = str2double(dateToken{1});
    else
        % If no date is found, set the value to 0
        fileDates(i) = 0;
    end
end

% Find the index of the most recent date
[~, idx] = max(fileDates);

% Get the most recent file
mostRecentFile = fullfile(directory, matchingFiles(idx).name);
end

% function mostRecentFile = findMostRecentFile(directory, token)
%     % Get a list of all files in the directory
%     files = dir(fullfile(directory, ['*' token '*']));
%
%     if isempty(files)
%         error('No files found containing the token "%s".', token);
%     end
%
%     % Extract the dates from the filenames
%     fileDates = zeros(length(files), 1);
%     for i = 1:length(files)
%         % Assume the date is at the end of the filename before the extension
%         [~, name, ~] = fileparts(files(i).name);
%
%         % Find the last set of 8 consecutive digits (assumed to be yyyymmdd)
%         dateToken = regexp(name, '\d{8}$', 'match');
%
%         if ~isempty(dateToken)
%             % Convert the dateToken string to a number for comparison
%             fileDates(i) = str2double(dateToken{1});
%         else
%             % If no date is found, set the value to 0
%             fileDates(i) = 0;
%         end
%     end
%
%     % Find the index of the most recent date
%     [~, idx] = max(fileDates);
%
%     % Get the most recent file
%     mostRecentFile = fullfile(directory, files(idx).name);
% end