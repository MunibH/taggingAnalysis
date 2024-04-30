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

% %% PHASE OF JAW
% 
% close all;
% clear jaw jawpos jawpos_norm
% 
% feats = {'jaw_ydisp_view1'};
% featmask = ismember(kin.featLeg,feats);
% 
% licktime = [-0.25 1]; % approximate first lick time
% ix = findTimeIX(obj.time,licktime);
% ix = ix(1):ix(2);
% % ix = 1:numel(obj.time);
% tm = obj.time(ix);
% 
% ct = 1;
% for trial = sort(randperm(obj.bp.Ntrials,1))
%     jawpos(:,ct) = squeeze( kin.dat(ix,trial,featmask) );
%     jawpos_norm(:,ct) = mySmooth(normalize(jawpos(:,ct),'range',[0 1]),7,'reflect');
%     closed_phase = 0; % phase of jaw when closed
%     [jaw.phase(:,ct), jaw.amp(:,ct)] = calculateJawPhase(jawpos(:,ct));
%     ct = ct + 1;
% end
% 
% % incix = findIncreaseIX(jawpos_norm);
% 
% % jaw.phase = mySmooth(jaw.phase,11,'reflect');
% % jaw.amp = mySmooth(jaw.amp,11,'reflect');
% 
% cm = colormap(flip(gray(numel(tm))));
% colorData = (jaw.phase - min(jaw.phase)) / (max(jaw.phase) - min(jaw.phase));
% close all
% 
% f = figure;
% f.Position = [139   410   560   420];
% ax1 = prettifyAxis(subplot(2,1,1));
% ax2 = prettifyAxis(subplot(2,1,2));
% % ax3 = prettifyAxis(subplot(3,1,3));
% % linkaxes([ax1,ax2,ax3],'x');
% hold(ax1,'on'); hold(ax2,'on'); %hold(ax3,'on')
% plot(ax1,obj.time(ix),jawpos,'k')
% s1=scatter(ax1,obj.time(ix),jawpos,20,colorData,'filled');
% colormap(cm)
% % plot(ax1,tm(incix), jawpos_norm(incix), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% plot(ax2,obj.time(ix),jaw.phase,'k');
% s2=scatter(ax2,obj.time(ix),jaw.phase,20,colorData,'filled');
% colormap(cm)
% % plot(ax2,tm(incix), jaw.phase(incix), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
% % plot(ax3,obj.time(ix),jaw.amp)
% xlabel(ax2,'Time from go cue (s)')
% ylabel(ax1,'Jaw Pos (px)')
% ylabel(ax2,'Phase (rad)')
% 
% 
% 
% 
% 
% gifname = 'test.gif';
% f = figure;
% f.Position = [720   410   560   420];
% ax = polaraxes;
% hold(ax,'on')
% % ax.RLim = [min(jawpos) max(jawpos)];
% ax.RLim = [0 1];
% ax.Title.String = ['Trial ' num2str(trial)];
% % cm = colormap(flip(gray(size(jaw.phase,1))));
% cm = repmat(s2.CData,1,3);
% for j = 1:numel(ix)
%     % ax.RLim = [min(jawpos) max(jawpos)];
%     ax.RLim = [0 1];
%     ps = polarscatter(jaw.phase(j,:),jaw.amp(j,:), [], cm(j,:), 'filled');
%     % ps = polarscatter(jaw.phase(j,:),jawpos(j,:), [], cm(j,:), 'filled');
%     for i = 1:numel(ps)
%         ps(i).SizeData = 35;
%     end
%     ax.ThetaTickLabel = arrayfun(@(x) sprintf('%.2f', x), deg2rad(ax.ThetaTick), 'UniformOutput', false);
% 
%     % Capture the frame
%     frame = getframe(f);
% 
%     % Convert the frame to an indexed image
%     im = frame2im(frame);
%     [indIm, cmap] = rgb2ind(im, 256);
% 
%     % Write the frame to the GIF file
%     if j == 1
%         % Write the first frame with the loop count and delay time
%         imwrite(indIm, cmap, gifname, 'gif', 'LoopCount', inf, 'DelayTime', 0.05);
%     else
%         % Append subsequent frames to the GIF file
%         imwrite(indIm, cmap, gifname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
%     end
% 
% 
%     drawnow;
% end
% 
% 
% 








