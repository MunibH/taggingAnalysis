function probeStr = generateTrialdatString(numProbes)
    % Initialize an empty array to store the strings
    probeStrings = strings(1, numProbes);
    
    % Loop over each probe and generate the string 'probeparams{i}.cluid'
    for i = 1:numProbes
        probeStrings(i) = sprintf('probeobj{%d}.trialdat ', i);
    end
    
    % Concatenate all the strings and enclose them in double quotes
    probeStr = sprintf('''{%s}''', strjoin(probeStrings));
end