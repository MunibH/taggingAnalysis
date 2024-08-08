clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%%
datapth = fullfile(utilspth,'tagVidData');
contents = dir(datapth);

v = VideoReader(fullfile(datapth, 'JPV13_2023-10-03_tag1_cam_0_trial3.mp4'));
load(fullfile(datapth,"data_structure_JPV13_2023-10-03.mat"))
load(fullfile(datapth,"SpkData_chans53-57.mat"))
load(fullfile(datapth,"SpkData_chans53-57_SortState-13-Oct-2023_done.mat"))

%%
close all


p.trial = 3;

p.vidtime = (1:v.NumFrames)./400;

p.window = 500 / 1000; % (ms) show this much data at a time
p.moveWin = 0.0025; % (ms) move forward in time this much per frame

p.in = clusters & spk.trial==p.trial;
p.out = ~clusters & spk.trial==p.trial;

p.stim = obj.sglx.laserTrigIX{p.trial} / obj.sglx.fs - 0.55;



f = figure; 
f.Position = [493   145   933   770];
f.Color = [0,0,0];
set(f,'InvertHardcopy','off')

ax1 = axes('Parent',f,... % spikes
    'Units','normalized',...
    'Position',[0 0.05 1 0.3],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);

% f1 = figure; % spikes
% f1.Position = [486    96   932   265];
% f1.Color = [0,0,0];
% set(f1,'InvertHardcopy','off')
% ax1 = gca;
ax1.Color = [0 0 0];
ax1.Visible = 'off';

ax2 = axes('Parent',f,... % video
    'Units','normalized',...
    'Position',[0.25 0.38 0.5 0.5],...
    'NextPlot','add',...
    'Visible','off',...
    'XTick',[],...
    'YTick',[]);

% f2 = figure; % video
% f2.Color = [0,0,0];
% set(f2,'InvertHardcopy','off')
% ax2 = gca;
ax2.Color = [0 0 0];
ax2.Visible = 'off';

start = 11;
stop = start + p.window;
% last = p.vidtime(end);
last = 15;

k = 1;
v.CurrentTime = 1;
while stop < last

    if k==1
        v.CurrentTime = start;
    end

    cla(ax1)
    cla(ax2)
    hold(ax1,'on');
    hold(ax2,'on');
    axtime = linspace(start,stop,1000);
    yline(ax1,0,'y-','LineWidth',1)
    line(ax1,[stop-0.1 stop-0.1],[50 100],'Color','w','LineWidth',2)
    text(ax1,stop-0.095,75,'200 mV','color','w','fontsize',20)
    xlim(ax1,[axtime(1),axtime(end)])
    ylim(ax1,[-100,100])

    % text time
    if mod(start,0.25)
        tstring = [num2str(round(start,2)) ' sec'];
        tt = text(ax1,axtime(end-200),-95,tstring);
        tt.Color = 'w';
        tt.FontSize = 20;
    end

    % plot in data that occurred during start:stop
    mask = (spk.trialtime+0.2)>=start & (spk.trialtime-0.2)<=stop;
    mask = mask & p.in;
    spktime = spk.trialtime(mask);
    spktrial = spk.trial(mask);
    spkwav = spk.wav(mask,:,1);
    for i = 1:numel(spktime)
        stm = spktime(i);
        swv = spkwav(i,:);
        xx = linspace(stm-0.05,stm+0.05,numel(swv));
        plot(ax1,xx,swv,'Color',[255, 184, 253]/255,'linewidth',2)
    end

    % plot out data that occurred during start:stop
    mask = (spk.trialtime+0.2)>=start & (spk.trialtime-0.2)<=stop;
    mask = mask & p.out;
    spktime = spk.trialtime(mask);
    spktrial = spk.trial(mask);
    spkwav = spk.wav(mask,:,1);
    for i = 1:numel(spktime)
        stm = spktime(i);
        swv = spkwav(i,:);
        xx = linspace(stm-0.03,stm+0.03,numel(swv));
        plot(ax1,xx,swv,'Color',[199, 199, 199]/255,'linewidth',2)
    end

    % plot stim
    mask = p.stim>=start & p.stim<=stop;
    if any(mask)
        stims = p.stim(mask);
        for i = 1:numel(stims)
            pch = patch(ax1,[stims(i),stims(i)+0.005,stims(i)+0.005,stims(i)],[-100 -100 100 100],[138, 239, 255]/255);
            pch.EdgeColor = 'none';
            pch.FaceAlpha = 0.6;
        end
    end

    % ensure axes look good
    xlim(ax1,[axtime(1),axtime(end)])
    ylim(ax1,[-100,100])

    xlabel(ax1,'Time (s)','fontsize',13)
    ylabel(ax1,'Volts','fontsize',13)



    % % plot video
    % Read video frames until available
    vidFrame = readFrame(v);
    ii = image(vidFrame, 'Parent', ax2);
    ax2.Visible = 'off';
    ax2.YDir = 'reverse';


    % update start and stop
    start = start + p.moveWin;
    stop = stop + p.moveWin;

    % draw and save
    drawnow
    s(k) = getframe(f);
    % s1(k) = getframe(ax1);
    % s2(k) = getframe(ax2);
    k = k + 1;
end


%% SAVE

outpth = fullfile(utilspth,'tagVidData');
vOut = VideoWriter(fullfile(outpth,'JPV13_2023-10-03_trial3.mp4'),'MPEG-4');
vOut.FrameRate = 35;
open(vOut)
for k = 1:numel(s)
    writeVideo(vOut,s(k))
end
close(vOut)

































