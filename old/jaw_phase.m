clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here

params.region = 'any'; % 'alm','tjm1','mc', 'any'
params.probeType = 'any'; % 'h2','np2','np1', 'any'

params.behav_only = 1;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
% meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
meta = loadJPV13(meta,datapth);

% subset meta
meta = meta(2);
% meta = subsetMetaByParams(meta,params);

%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

me = loadMotionEnergy(obj, meta, params, datapth);

kin = getKinematics(obj, me, params);

% TAGGED UNIT META

tag.nTag = numel(obj.tag);
tag.cluid.clu = [obj.tag(:).cluid]; % where tagged units are in obj.clu
tag.cluid.obj = find(ismember(params.cluid,tag.cluid.clu))'; % where in params.cluid, trialdat, psth

%% OBJVIS
clc,close all
addpath(genpath(fullfile(utilspth,'objvis')));
objvis(meta,params,obj,tag,kin)

%%

feats = {'jaw_phase'};
featmask = ismember(kin.featLeg,feats);


tix = findTimeIX(obj.time,[-0.2 1],1);
jp = kin.dat(tix,:,featmask);

f = figure;
ax = prettifyAxis(gca);
imagesc(obj.time(tix),1:obj.bp.Ntrials,jp')

Fs = 1/params.dt;
T = 1/Fs;
L = numel(tix);
t = (0:L-1)*T;

Y = fft(jp(:,1)); % fourier coef for a single trial

f = figure;
ax = prettifyAxis(gca);
plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y)),"LineWidth",3)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f (Hz)")
ylabel("|fft(X)|")



%% PHASE OF JAW

close all;
clear jaw jawpos jawpos_norm

% trial = 113;
trial = 139;
% trial = randperm(obj.bp.Ntrials,1);

feats = {'jaw_ydisp_view1'};
featmask = ismember(kin.featLeg,feats);

frametimes = obj.traj{1}(trial).frameTimes - obj.bp.ev.goCue(trial);

times = [-0.1 1.5];
% times = [0 0.3]; 
ix = findTimeIX(frametimes,times,1);
% ix = 1:numel(obj.time);
tm = frametimes(ix);

fs = 400;

% [b,a] = butter(3,[0.1/fs/2 40/fs/2], 'bandpass'); % [0.1/fs/2 40/fs/2] with islocalminima on phase works ok at finding jaw closing times
[b,a] = butter(3,50/fs/2, 'low'); % [0.1/fs/2 40/fs/2] with islocalminima on phase works well at finding jaw closing times


jawpos = obj.traj{1}(trial).ts(ix,2,4); % jaw y pos in pixels
% jawpos = mySmooth(jawpos,31,'reflect');
jawpos = filtfilt(b,a,jawpos); % bandpass filter, want narrow band signal to avoid sharp transitions in phase extracted by hilbert transform

% jawpos = squeeze( kin.dat(ix,trial,featmask) );
[jaw.phase, jaw.amp] = calculateJawPhase(jawpos);

% incix = findIncreaseIX(jawpos_norm);

% jaw.phase = mySmooth(jaw.phase,11,'reflect');
% jaw.amp = mySmooth(jaw.amp,11,'reflect');

cm1 = colormap(flip(copper(ceil(numel(tm)/2))));
cm2 = flip(cm1);
cm = cat(1,cm1,cm2);
colorData = (jaw.phase - min(jaw.phase)) / (max(jaw.phase) - min(jaw.phase));
close all

% find local minima in phase to get individual licks
TF = islocalmin(jaw.phase);

f = figure;
f.Position = [139   410   560   420];
ax1 = prettifyAxis(subplot(2,1,1));
ax2 = prettifyAxis(subplot(2,1,2));
% ax3 = prettifyAxis(subplot(3,1,3));
% ax3 = prettifyAxis(subplot(3,1,3));
% linkaxes([ax1,ax2,ax3],'x');
hold(ax1,'on'); hold(ax2,'on'); %hold(ax3,'on')
plot(ax1,tm,jawpos,'k')
s1=scatter(ax1,tm,jawpos,20,colorData,'filled');
% s11=scatter(ax1,tm(TF),jawpos(TF),40,'red','filled');
colormap(cm)
% plot(ax1,tm(incix), jawpos_norm(incix), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
plot(ax2,tm,jaw.phase,'k');
s2=scatter(ax2,tm,jaw.phase,20,colorData,'filled');
% s22=scatter(ax2,tm(TF),jaw.phase(TF),40,'red','filled');
colormap(cm)
% plot(ax2,tm(incix), jaw.phase(incix), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% plot(ax3,obj.time(ix),jaw.amp)

xlabel(ax2,'Time from go cue (s)')
ylabel(ax1,'Jaw Pos (px)')
ylabel(ax2,'Phase (rad)')
% jaw.phasedt = abs(gradient(jaw.phase));
% [pks,locs] = findpeaks(jaw.phasedt,'MinPeakHeight',0.05,'MinPeakDistance',30);
% plot(ax3,tm, jaw.phasedt,'k');
% s3=scatter(ax3,tm,jaw.phasedt,20,colorData,'filled');
% s4 = scatter(ax3,tm(locs),jaw.phasedt(locs),20,'red','filled');
% colormap(cm)




%%
gifname = 'test.gif';
f = figure;
f.Position = [720   410   560   420];
ax = polaraxes;
hold(ax,'on')
ax.RLim = [min(jawpos) max(jawpos)];
% ax.RLim = [0 1];
ax.Title.String = ['Trial ' num2str(trial)];
% cm = colormap(flip(gray(size(jaw.phase,1))));
% cm = repmat(s22.CData,1,3);
for j = 400:600%1:numel(ix)
    ax.RLim = [min(jawpos) max(jawpos)];
    % ax.RLim = [0 1];
    % ps = polarscatter(jaw.phase(j,:),jaw.amp(j,:), [], cm(j,:), 'filled');
    ps = polarscatter(jaw.phase(j,:),jawpos(j,:), [], cm(j,:), 'filled');
    for i = 1:numel(ps)
        ps(i).SizeData = 35;
    end
    ax.ThetaTickLabel = arrayfun(@(x) sprintf('%.2f', x), deg2rad(ax.ThetaTick), 'UniformOutput', false);

    % Capture the frame
    frame = getframe(f);

    % Convert the frame to an indexed image
    im = frame2im(frame);
    [indIm, cmap] = rgb2ind(im, 256);

    % Write the frame to the GIF file
    if j == 1
        % Write the first frame with the loop count and delay time
        imwrite(indIm, cmap, gifname, 'gif', 'LoopCount', inf, 'DelayTime', 0.05);
    else
        % Append subsequent frames to the GIF file
        imwrite(indIm, cmap, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end


    drawnow;
end











