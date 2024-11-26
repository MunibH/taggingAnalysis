clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'rastermap')));
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

params.smooth = 11;

params.dt = 1/50;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
rastermappth = fullfile(datapth,'rastermap');
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)
meta = meta(1);

%% save rastermap data

for isess = 1:numel(meta)

    thismeta = meta(isess);

    [obj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(obj,sesspar,thismeta);

    %
    extractRasterMapData(rastermappth,thismeta,obj,params,sesspar,tag)

    break

end



%% Helper functions

function extractRasterMapData(rastermappth,thismeta,obj,params,sesspar,tag)

%%
sessionEndTime = obj.sglx(1).fileStart(end) + (obj.sglx(1).Nsamp(end)./obj.sglx(1).fs); % session end in seconds
edges = -0.5:params.dt:sessionEndTime;


% offset = (obj.sglx.fileStart - obj.sglx.bitFileOffset);
offset = obj.sglx.fileStart;
gocue = obj.bp.ev.goCue + offset;
trialstart = obj.bp.ev.bitStart + offset;
sample = obj.bp.ev.sample + offset;
delay = obj.bp.ev.delay + offset;
dt = params.dt;
tedges = edges;

nProbes = numel(sesspar.cluid);
nClu = sum(cell2mat(cellfun(@numel, sesspar.cluid,'uni',0)));
spks = nan(nClu,numel(edges));
tagged = false(nClu,1);

ct = 1;
for iprobe = 1:nProbes
    for iclu = 1:numel(sesspar.cluid{iprobe})
        thisclu = sesspar.cluid{iprobe}(iclu);

        N = histc(obj.clu{iprobe}(thisclu).tm, edges);
        N = mySmooth(N./params.dt,params.smooth,params.bctype);
        % spks(ct,:) = zscore(N);
        spks(ct,:) = N; % rastermap zscores already
        if strcmpi(obj.clu{iprobe}(thisclu).quality,'tagged')
            tagged(ct) = true;
        end
        ct = ct + 1;
    end
end

psth = cat(2,obj.psth{1}(:,:,[2 3]),obj.psth{2}(:,:,[2,3]));

anm = thismeta.anm;
sessiondate = thismeta.date;
time = obj.time;

savepth = rastermappth;
savefn = [thismeta.anm '_' thismeta.date '_RasterMapData.mat'];
if ~exist(savepth)
    mkdir(savepth)
end
save(fullfile(savepth,savefn),'anm','sessiondate','spks','gocue','trialstart','sample','delay','dt','tedges','tagged','psth','time')


disp(['Saved: ' fullfile(savepth,savefn)])





end


%%
function events = getFirstHighBitcodeIX(meta)
Ntrials = numel(meta.sglxfns);
events = NaN.*zeros(Ntrials, 1);

firstIX = 1;
noUse = true;
while noUse
    if isempty(meta.bitcode.highIX{firstIX})
        firstIX = firstIX+1;
    else
        noUse = false;
        events(firstIX) = (meta.bitcode.highIX{firstIX}(1))./meta.fs;
    end
end

for i = firstIX:Ntrials-1
    if ~isempty(meta.bitcode.highIX{i+1})
        % number of samples recorded up to previous trial plus index of current trial go cue
        events(i+1) = (sum(meta.Nsamp(1:i)) + meta.bitcode.highIX{i+1}(1)) ./ meta.fs;
    end
end
events = events(~isnan(events));

end % getGoCue

