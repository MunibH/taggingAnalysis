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
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
% params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};
params.qm.quality = {'tagged'};

% params.alignEvent = 'lastLick';
% params.tmin = -4;
% params.tmax = 2;

params.behav_only = 0;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% subset meta (TODO)

% meta = subsetMetaByParams(meta,params);


%% LOAD DATA

for isess = 1:numel(meta)
    disp(' ')
    disp(['Session ' num2str(isess) '/' num2str(numel(meta))])
    disp(' ')
    [sessobj,sesspar] = loadSessionData(meta(isess),params);

    tag(isess) = getTagFromObj(sessobj,sesspar,meta(isess));
    break
end

%% SAVE TAG

fpth = 'C:\Users\munib\Documents\Economo-Lab\data\tagged';
% fn = 'AllTagged_GoCue_20240817.mat';
fn = 'AllTagged_LastLick_20240817.mat';
% fn = 'AllTagged_FirstLick_20240817.mat';
% save(fullfile(fpth,fn),'tag','params','-v7.3')







