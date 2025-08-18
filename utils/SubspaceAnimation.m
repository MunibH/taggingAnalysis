% anti_correlated_poisson_demo.m
% Simulate two neurons with perfectly anti‑correlated firing rates,
% generate inhomogeneous Poisson spikes, and animate rasters + FR state‑space.

clear; close all;

%% Simulation parameters
T    = 10;        % total duration (s)
dt   = 0.05;     % time step (s)
t    = 0:dt:T;    % time vector
base = 20;        % baseline firing rate (Hz)
amp  = 15;        % modulation amplitude (Hz)
f    = 0.5;       % modulation frequency (Hz)
mx = 50;
mn = 0;
base = 25;

n = length(t)/4;

%% 1) Generate time‑varying firing rates
% Neuron 1: sinusoidal around 'base'

rates1 = cat(2, linspace(base,mx,n), linspace(mx,base,n), linspace(base,mn,n), linspace(mn,base,n+1) );

rates2 = cat(2, linspace(base,mx,n), linspace(mx,base,n), linspace(base,mn,n), linspace(mn,base,n+1) );

%% 2) Simulate inhomogeneous Poisson spike trains
% At each dt, spike probability = rate*dt
spk1 = (rand(size(t)) < abs(rates1) * dt) & (rand(size(t)) < 0.75);
spk2 = (rand(size(t)) < abs(rates2) * dt) & (rand(size(t)) < 0.75);

%% 3) Estimate firing rates via sliding‑window smoothing (100 ms)
% win    = round(0.1/dt);
% fr1    = conv(spk1, ones(1,win)/win, 'same') / dt;
% fr2    = conv(spk2, ones(1,win)/win, 'same') / dt;

fr1 = rates1;
fr2 = rates2;

%% 4) Set up figure
figure('Color','w','Position',[100 100 600 800]);

% Top: spike rasters
ax1 = subplot(2,1,1);
h1  = plot(ax1, nan, 1, '.k','MarkerSize',10); hold on;
h2  = plot(ax1, nan, 2, '.b','MarkerSize',10);
ylim([0.5 2.5]);
yticks([1 2]); yticklabels({'Neuron 1','Neuron 2'});
xlabel('Time (s)'); title('Spike Rasters');
xlim([0 T]);

% Bottom: FR state‑space
ax2 = subplot(2,1,2);
hTraj   = plot(ax2, nan, nan, 'LineWidth',1.5); hold on;
hPoint  = plot(ax2, nan, nan, 'ro','MarkerFaceColor','r');
% Null subspace line: FR2 = 2*base – FR1
lims = [mn mx];
axis equal; xlim(lims); ylim(lims);
xlabel('firing rate neuron 1'); ylabel('firing rate neuron 2');


%% 5) Animate
for i = 1:length(t)
    % Update rasters up to current time
    ti = t(i);
    set(h1, 'XData', t(spk1 & (t<=ti)), 'YData', ones(1,sum(spk1 & (t<=ti))));
    set(h2, 'XData', t(spk2 & (t<=ti)), 'YData', ones(1,sum(spk2 & (t<=ti)))*2);


    
    % Update FR trajectory and current point
    set(hTraj, 'XData', fr1(1:i), 'YData', fr2(1:i));
    set(hPoint,'XData', fr1(i),   'YData', fr2(i));
    
    drawnow;
    % pause(0.01);
end
