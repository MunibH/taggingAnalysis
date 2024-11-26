clear,clc,close all

% Parameters
x = linspace(-5, 5, 100)'; % Input space
sigma = 1;              % RBF kernel width
len = 0.01;                  % Kernel scale

% Compute the RBF kernel matrix
K = rbf(x, sigma, len);

% Mean and covariance for the GP
mu = zeros(size(x));      % GP mean
L = chol(K + 1e-6 * eye(size(K))); % Add jitter for numerical stability

% Initialize figure
f = figure;
ax = gca;
hold on;
grid on;
xlabel('Input x');
ylabel('f(x)');
title('Samples from a Gaussian Process with RBF Kernel');
set(ax, 'Color', [0.8 0.8 0.8]); % Light background for fading effect

% Draw samples continuously
prev_x = [];
prev_y = [];
while true
    % Draw a random sample from the GP
    f = mu + L' * randn(size(x));
    
    if ~isempty(prev_x)
        cla(gca);
        plot(prev_x,prev_y,'LineWidth',1.3,'Color',[0.3,0.3,0.3]);
    end

    % Plot the new sample
    h = plot(x, f, 'LineWidth', 2, 'Color', 'b');
    
    % Pause and allow interaction
    % pause(0.05);
    drawnow;

    prev_x = cat(2,prev_x,h.XData');
    prev_y = cat(2,prev_y,h.YData');
end

% RBF kernel function
function K = rbf(x, sigma, len)
    x = x(:);
    dists = pdist2(x, x).^2;
    K = len^2 * exp(-dists / (2 * sigma^2));
end
