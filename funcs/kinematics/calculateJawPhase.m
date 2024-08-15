function [phase,amplitude,analyticalSignal,jawpos_filt] = calculateJawPhase(jawpos,Fs)
%%
% % first want to filter out unwanted frequencies
% % Parameters
% Fs = Fs;  % Sampling frequency (in Hz, replace with your actual sampling rate)
% Fc = 15;    % Cutoff frequency (in Hz), anything above this is ousted
% 
% % Normalized cutoff frequency (between 0 and 1)
% Wn = Fc / (Fs / 2);
% 
% % Design a 4th-order Butterworth low-pass filter
% [b, a] = butter(4, Wn, 'low');
% % filter
% % jaw = mySmooth(jawpos,21,'reflect');
% % Fill missing values using a moving average

%%
% % Y = fft(jawpos);
% % f = (0:length(Y)-1) * (400 / length(Y));
% % plot(f, abs(Y),'LineWidth',2);
% % xlabel('Frequency (Hz)');
% % ylabel('Magnitude');
% % frequency content of jaw has significant peaks from 1-15 Hz
% % let's bandpass 

% Parameters
Fs = Fs;  % Sampling frequency (in Hz, replace with your actual sampling rate)
Fc1 = 0.1;    % Lower cutoff frequency (in Hz)
Fc2 = 15;   % Upper cutoff frequency (in Hz)

% Normalized cutoff frequencies (between 0 and 1)
Wn = [Fc1 Fc2] / (Fs / 2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(4, Wn, 'bandpass');


x_filled = fillmissing(jawpos, 'movmean', 10);  % Replace window_size with an appropriate value
jawpos_filt = filtfilt(b,a,x_filled);

% figure; hold on; 
% yyaxis left
% plot(jawpos,'LineWidth',2); 
% ylabel('jaw pos')
% yyaxis right
% plot(jawpos_filt,'LineWidth',2)
% ylabel('filtered')

% Compute the analytical signal using the Hilbert transform
analyticalSignal = hilbert(jawpos_filt);
% Extract the phase from the analytical signal
phase = angle(analyticalSignal);
% amplitude envelope
% amplitude = normalize(jawpos,'range',[0 1]);
amplitude = abs(analyticalSignal);
% % % 
% % % % Define a threshold for oscillation detection
% % % threshold = some_value;  % Set this based on your signal characteristics
% % % 
% % % % Ignore phase where amplitude is below threshold
% % % phase(amplitude_envelope < threshold) = NaN;  % Optional step

%%

end