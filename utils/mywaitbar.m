function mywaitbar(i, totalIterations)

if i == 1
    % Initialize the text-based loading bar
    fprintf('Progress: 0%% [');
end
for i = 1:totalIterations
    % Update the text-based loading bar
    progress = i / totalIterations;
    barLength = 50;  % Length of the loading bar
    barString = repmat('=', 1, round(progress * barLength));
    emptyString = repmat(' ', 1, barLength - length(barString));
    fprintf('\b\b\b\b\b\b\b\b\b\b');  % Clear the previous loading bar
    fprintf('%s%s] %d%%', barString, emptyString, round(progress * 100));
end
if i==totalIterations
    fprintf('\n');  % Start a new line after the loading bar is complete
end

end