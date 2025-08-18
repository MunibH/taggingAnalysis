% three_neuron_state_space_with_line.m
% Simulate three neurons with a smooth 3D firing‑rate trajectory,
% generate inhomogeneous Poisson spikes, and animate:
%  1) Spike rasters in the top subplot
%  2) 3D neural trajectory with colored points and an evolving line in the bottom subplot

clear; close all;

rng(25)

%% Simulation parameters
T    = 4;        % total time (s)
dt   = 0.05;      % time step (s)
t    = 0:dt:T;    % time vector
f    = 0.2;       % oscillation frequency (Hz)
base = 20;        % baseline rate (Hz)
amp  = 15;        % modulation amplitude (Hz)

%% Generate smooth 3D firing‑rate trajectories (phase‐shifted sinusoids)
rates1 = 20 + amp * sin(2*pi*f*t);
rates2 = 10 + amp * sin(2*pi*f*t + 2*pi/3);
% rates3 = 50 + amp * -cos(2*pi*f*t + 4*pi/3);
rates3 = linspace(500,100,length(t));

%% Simulate inhomogeneous Poisson spikes
spk1 = rand(size(t)) < rates1 * dt  & (rand(size(t)) < 0.75);
spk2 = rand(size(t)) < rates2 * dt  & (rand(size(t)) < 0.75);
spk3 = rand(size(t)) < rates3 * dt  & (rand(size(t)) < 0.75);

%% Smooth spikes back into rate estimates (100 ms window)
win = gausswin(61);
fr1 = conv(spk1, win, 'same') / dt;
fr2 = conv(spk2, win, 'same') / dt;
fr3 = conv(spk3, win, 'same') / dt;

%% Colormap for time
cmap = plasma(length(t));
cmapproj = gray(length(t));

spikecmap = linspecer(3);

%% Set up figure
figure('Color','w','Position',[7 374 1484 481]);

% Top subplot: spike rasters
ax1 = subplot(1,2,1);
ax1.Position = [0.1 0.3 0.4 0.3];
ax1.Box = 'off';
% h1 = plot(ax1, NaN,1,'.k','MarkerSize',10); hold on;
% h2 = plot(ax1, NaN,2,'.b','MarkerSize',10);
% h3 = plot(ax1, NaN,3,'.r','MarkerSize',10);
h1 = plot(ax1, NaN,1,'|','MarkerSize',30,'LineWidth',2,'Color',spikecmap(1,:)); hold on;
h2 = plot(ax1, NaN,2,'|','MarkerSize',30,'LineWidth',2,'Color',spikecmap(2,:));
h3 = plot(ax1, NaN,3,'|','MarkerSize',30,'LineWidth',2,'Color',spikecmap(3,:));
ylim([0.7 3.3]); yticks(1:3); 
yticklabels({'Neuron 1','Neuron 2','Neuron 3'});
ax1.TickDir = 'none';
ax1.YAxis.LineWidth = 1;
ax1.YAxisLocation = 'left';
xlabel('Time'); 
xlim([0 T]);
ax1.XTick = [];
ax1.FontSize = 20;
ax1.XColor = 'none';
xlabel('Time','Color','k'); 

% Bottom subplot: 3D state‑space
ax2 = subplot(1,2,2);
hLine  = plot3(ax2, NaN,NaN,NaN, '-', 'LineWidth',1.5, 'Color', [0.3 0.3 0.3]); hold on;
hTrail = scatter3(ax2, NaN,NaN,NaN, 100, NaN, 'filled');
hPoint = scatter3(ax2, NaN,NaN,NaN,300, NaN, 'filled');

hLineProj  = plot3(ax2, NaN,NaN,NaN, '-', 'LineWidth',3, 'Color', [0.3 0.3 0.3]); hold on;
hTrailProj = scatter3(ax2, NaN,NaN,NaN, 100, NaN, 'filled');
hPointProj = scatter3(ax2, NaN,NaN,NaN,300, NaN, 'filled');

buf = 50;
axis(ax2, [min(fr1)-buf max(fr1)+buf...
    min(fr2)-buf max(fr2)+buf...
    min(fr3)-buf max(fr3)+buf]);
view(107.6048,22.25);
% ax2.XColor = 'none';
% ax2.YColor = 'none';
% ax2.ZColor = 'none';
xlabel({'firing rate' ,'neuron 1'}); 
ylabel('neuron 2'); 
zlabel('neuron 3');
ax2.FontSize = 20;
ax2.XTickLabel = [];
ax2.YTickLabel = [];
ax2.ZTickLabel = [];
grid(ax2,'on')
ax2.LineWidth = 1;


%% Animation
for i = 1:length(t)
    ti = t(i);
    
    % Update spike rasters
    set(h1, 'XData', t(spk1 & t<=ti), 'YData', ones(1,sum(spk1 & t<=ti)));
    set(h2, 'XData', t(spk2 & t<=ti), 'YData', ones(1,sum(spk2 & t<=ti))*2);
    set(h3, 'XData', t(spk3 & t<=ti), 'YData', ones(1,sum(spk3 & t<=ti))*3);
    
    % Indices up to current time
    idx = 1:i;
    
    % Update evolving line
    set(hLine, 'XData', fr1(idx), 'YData', fr2(idx), 'ZData', fr3(idx));
    
    % Update colored trail points
    set(hTrail, 'XData', fr1(idx), 'YData', fr2(idx), 'ZData', fr3(idx), ...
                'CData', cmap(idx,:));
    
    % Update current point
    set(hPoint, 'XData', fr1(i), 'YData', fr2(i), 'ZData', fr3(i), ...
                'CData', cmapproj(i,:));

    % proj
    % set(hLineProj, 'XData', fr1(idx), 'YData', zeros(size(fr1(idx))), 'ZData', fr3(idx)); 
    % set(hTrailProj, 'XData', fr1(idx), 'YData', zeros(size(fr1(idx))), 'ZData', fr3(idx), ...
    %             'CData', cmapproj(idx,:));
    % set(hPointProj, 'XData', fr1(i), 'YData', zeros(size(fr1(i))), 'ZData', fr3(i), ...
    %             'CData', cmapproj(i,:));

    set(hLineProj, 'XData', fr1(idx), 'YData', fr2(idx), 'ZData', min(fr3)*ones(size(fr1(idx)))-buf); 
    set(hTrailProj, 'XData',fr1(idx), 'YData', fr2(idx), 'ZData', min(fr3)*ones(size(fr1(idx))-buf), ...
                'CData', cmapproj(idx,:));
    set(hPointProj, 'XData', fr1(i), 'YData', fr2(i), 'ZData', min(fr3)-buf, ...
                'CData', [0 0 0]);
    
    drawnow;
    pause(0.075);
end
