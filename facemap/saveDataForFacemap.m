clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'facemap')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

params.behav_only = 0;

params.smooth = 51;

params.dt = 1/75;


%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)


%% LOAD DATA

thismeta = meta(1);

[obj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(obj,sesspar,thismeta);

%% load motion svd facemap results

% load('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data\MAH24_2024-06-11_proc_v2.mat');
load('C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data\JPV11_2023-06-16_proc.mat');
clearvars -except motSVD_0 movSVD_0 obj sesspar meta thismeta tag params datapth 

%% get number of frames in behav trials

clear nFrames frameTimes frameTrialStartTimes frameTrialEndTimes nCumFrames
for i = 1:obj.bp.Ntrials
    frameTimes{i} = obj.traj{1}(i).frameTimes - 0.5;
    nFrames(i) = numel(frameTimes{i});
    frameTrialStartTimes(i) = frameTimes{i}(1);
    frameTrialEndTimes(i) = frameTimes{i}(end);
end
nCumFrames = cumsum(nFrames);

%% temporally align motion svds and neural data

clear neuralData motSVD

sv2use = 1:100; % for each of motion and movie
neuralData = cell(1);
behavData = cell(1);
cluct = 1;
for iprobe = 1:numel(sesspar.probe)
    clus = sesspar.cluid{iprobe};
    for iclu = 1:numel(clus)
        thisclu = obj.clu{iprobe}(clus(iclu));
        c = 1;
        for itrial = 1:obj.bp.Ntrials
            onTrial = thisclu.trial==itrial;

            s = frameTrialStartTimes(itrial);
            e = frameTrialEndTimes(itrial);
            edges = s:params.dt:e;
            tvec = edges + params.dt/2;
            tvec = tvec(1:end-1);

            N = histc(thisclu.trialtm(onTrial), edges);
            neuralData{c}(:,cluct) = mySmooth(N(1:end-1),params.smooth,params.bctype);
            
            % resample motion sv's according to neural time
            ft = frameTimes{itrial};
            if iprobe==1 && iclu==1
                if itrial==1
                    iframes = 1:nCumFrames(itrial);
                else
                    iframes = nCumFrames(itrial-1)+1:nCumFrames(itrial);
                end
                motSVD = interp1(frameTimes{itrial},motSVD_0(iframes,sv2use),tvec); % interp1(old_time,me,new_time);
                movSVD = interp1(frameTimes{itrial},movSVD_0(iframes,sv2use),tvec); % interp1(old_time,me,new_time);
                behavData{c} = cat(2,motSVD,movSVD);
            end

            c = c + 1;
        end
        cluct = cluct + 1;
    end
end

%% vizualize if neural and kin data are aligned
% 
% close all
% f = figure;
% ax = gca;
% hold on;
% 
% cluid = 71; % probe 1 tagged is 616:621
% for itrial = 1:obj.bp.Ntrials
%     cla(ax);
% 
%     thisN = neuralData{itrial}(:,cluid);
%     thisM = motSVD{itrial}(:,1);
% 
%     edges = frameTrialStartTimes(itrial):params.dt:frameTrialEndTimes(itrial);
%     tvec = edges + params.dt/2;
%     tvec = tvec(1:end-1);
% 
%     thisN = normalize(thisN,'range',[0,1]);
%     thisM = normalize(thisM,'range',[0,1]);
% 
%     plot(tvec,thisN);
%     plot(tvec,thisM);
% 
%     tstart = obj.bp.ev.bitStart(itrial);
%     sample = obj.bp.ev.sample(itrial);
%     delay = obj.bp.ev.delay(itrial);
%     gocue = obj.bp.ev.goCue(itrial);
% 
%     plot([tstart,tstart],ax.YLim,'k--')
%     plot([sample,sample],ax.YLim,'k--')
%     plot([delay,delay],ax.YLim,'k--')
%     plot([gocue,gocue],ax.YLim,'k--')
% 
%     title([num2str(itrial) ' ' num2str(obj.bp.hit(itrial))])
%     pause
% end


%% save motion SV and neural data to .npy

addpath(genpath('C:\npy-matlab\npy-matlab'))

% nTimeEachTrial
nTimeEachTrial = cell2mat(cellfun(@(x) size(x,1), neuralData, 'UniformOutput',false));

% motion SV
videoSVD = cat(1,behavData{:}); % (time*trials,SVs)

% neural activity
neuralActivity = cat(1,neuralData{:}); % (time*trials,neurons)

dt = params.dt;

% tosave = {nTimeEachTrial,motionSVD,neuralActivity};

savepth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\facemap\data';
savefn = [thismeta.anm '_' thismeta.date '_FacemapData.mat'];
save(fullfile(savepth,savefn),'nTimeEachTrial','videoSVD','neuralActivity','dt')


disp('READY TO RUN FACEMAP-NN')




%% plot top 5 videoSVD

% close all
% 
% trialStarts = cumsum(nTimeEachTrial);
% 
% ix = 1:800;
% ix = 4000:6000;
% tstarts = trialStarts(trialStarts>=ix(1) & trialStarts<=ix(end)) .* params.dt;
% 
% 
% f = figure;
% f.Position = [372   462   655   448];
% ax = prettifyAxis(gca);
% cmap_sweep(10, mako); % set colororder for plot
% hold on;
% for i = 1:5
%     plot((ix).*params.dt, videoSVD(ix,i) + 100*(i-1) ,'LineWidth',0.1)
% end
% xlabel('Time (s)')
% ylabel('Top 5 ME PCs')
% xlim([ix(1).*params.dt ix(end)*params.dt])
% for i = 1:numel(tstarts)
%     plot([tstarts(i) tstarts(i)], ax.YLim,'-', 'Color',[0.3 0.3 0.3],'LineWidth',0.1)
% end
% 
% f = figure;
% f.Position = [372   462   655   448];
% ax = prettifyAxis(gca);
% cmap_sweep(10, rocket); % set colororder for plot
% hold on;
% ct = 1;
% for i = 101:106
%     plot((ix).*params.dt, videoSVD(ix,i) + 900*(ct-1) ,'LineWidth',0.1)
%     ct = ct + 1;
% end
% xlabel('Time (s)')
% ylabel('Top 5 Video PCs')
% xlim([ix(1).*params.dt ix(end)*params.dt])
% for i = 1:numel(tstarts)
%     plot([tstarts(i) tstarts(i)], ax.YLim,'-', 'Color',[0.3 0.3 0.3],'LineWidth',0.1)
% end



