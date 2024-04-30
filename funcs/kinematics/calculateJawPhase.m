function [phase,amplitude,analyticalSignal] = calculateJawPhase(jawpos)
    % jawpos is a 1d time series

    % Normalize the jaw movements
    jawpos = (jawpos - mean(jawpos)) / std(jawpos);

    % mean_x = mean(jawpos(:, 1));
    % mean_y = mean(jawpos(:, 2));
    % std_x = std(jawpos(:, 1));
    % std_y = std(jawpos(:, 2));
    % jawpos(:, 1) = (jawpos(:, 1) - mean_x) / std_x;
    % jawpos(:, 2) = (jawpos(:, 2) - mean_y) / std_y;


    % Compute the analytical signal using the Hilbert transform
    analyticalSignal = hilbert(jawpos);

    % amplitude = abs(analyticalSignal);
    % amplitude = sqrt(jawpos(:, 1).^2 + jawpos(:, 2).^2);

    % Calculate the angle (phase) from the jaw movements
    % phase = atan2(jawpos(:, 2), jawpos(:, 1));
    % Extract the phase from the analytical signal
    phase = angle(analyticalSignal);
    
    % Convert the angle from [-pi, pi] to [0, 2*pi]
    % phase = mod(phase,2*pi);
    % phase(phase < 0) = phase(phase < 0) + 2 * pi;
    
    % Adjust the phase if necessary (e.g., to set jaw closed at 0 and 2*pi)
    % Find the angle when the jaw is closed (hard coded and passed in as arg for now)
    
    % Subtract the closed angle from all phases to adjust
    % phase = mod(phase - closed_phase, 2 * pi);


    % % Detect peaks and troughs using findpeaks function
    % [peaks, peakLocs] = findpeaks(amplitude); % Peak values and their indices
    % [troughs, troughLocs] = findpeaks(-jawpos); % Trough values and their indices (negative values)
    % 
    % % Calculate the average phase at peaks and troughs
    % avgPhaseAtPeaks = mean(phase(peakLocs));
    % avgPhaseAtTroughs = mean(phase(troughLocs));
    % 
    % % Calculate the adjustment needed for peaks and troughs
    % % Adjust the average phase to be pi at peaks and 0/2pi at troughs
    % peakAdjustment = pi - avgPhaseAtPeaks
    % troughAdjustment = -avgPhaseAtTroughs; % since 0 or 2pi is the same in phase
    % 
    % phase = mod(phase + peakAdjustment, 2 * pi);

    amplitude = normalize(jawpos,'range',[0 1]);

end