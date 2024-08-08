function startIncreaseIdx = findIncreaseIX(periodic_signal)
% Calculate the derivative of the signal (approximation)
dx = diff(periodic_signal);

% Identify points where the derivative changes from negative to positive
% This indicates the start of an increase
% startIncreaseIdx = find(dx(1:end - 1) < 0 & dx(2:end) >= 0) + 1;

% Alternatively, you can use findpeaks with the negative of the signal to find local minima
[minValues, startIncreaseIdx] = findpeaks(-periodic_signal);

end