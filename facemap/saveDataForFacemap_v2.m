clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'facemap')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath('C:\npy-matlab\npy-matlab'))

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

params.nSVs = 100; % each of motion and movie (so 2x this number)

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
facemapsvdpth = fullfile(datapth,'Facemap');
meta = [];

meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)


%% save facemap data

for isess = 1:numel(meta)

    thismeta = meta(isess);

    [obj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(obj,sesspar,thismeta);

    % load facemap results, extract nSVs, align SVs and neural data, save
    extractFacemapData(facemapsvdpth,thismeta,obj,params,sesspar,tag)


end



%% Helper functions

function extractFacemapData(facemapsvdpth,thismeta,obj,params,sesspar,tag)
% load facemap svd results
contents = dir(facemapsvdpth);
files = {contents(:).name}';
mask = contains(files,thismeta.anm) & contains(files,thismeta.date) & contains(files,'.mat');
load(fullfile(facemapsvdpth,files{mask}))
clearvars -except motSVD_0 movSVD_0 obj sesspar meta thismeta tag params datapth meta facemapsvdpth

% get number of frames in behav trials
clear nFrames frameTimes frameTrialStartTimes frameTrialEndTimes nCumFrames
for i = 1:obj.bp.Ntrials
    frameTimes{i} = obj.traj{1}(i).frameTimes - 0.5;
    nFrames(i) = numel(frameTimes{i});
    frameTrialStartTimes(i) = frameTimes{i}(1);
    frameTrialEndTimes(i) = frameTimes{i}(end);
end
nCumFrames = cumsum(nFrames);


% temporally align motion svds and neural data

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

% save motion SV and neural data to .npy

% nTimeEachTrial
nTimeEachTrial = cell2mat(cellfun(@(x) size(x,1), neuralData, 'UniformOutput',false));

% motion SV
videoSVD = cat(1,behavData{:}); % (time*trials,SVs)

% neural activity
neuralActivity = cat(1,neuralData{:}); % (time*trials,neurons)

dt = params.dt;

% tosave = {nTimeEachTrial,motionSVD,neuralActivity};

savepth = facemapsvdpth;
savefn = [thismeta.anm '_' thismeta.date '_FacemapData.mat'];
save(fullfile(savepth,savefn),'nTimeEachTrial','videoSVD','neuralActivity','dt')


disp(['Saved: ' fullfile(savepth,savefn)])


end

