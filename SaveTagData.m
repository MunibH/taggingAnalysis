clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% process all tagged units and save to single .mat file

%% PARAMETERS

params = defaultParams();

% % specify changes here
params.alignEvent = 'goCue';
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
savepath = 'C:\Users\munib\Documents\Economo-Lab\data\curated\PTNs';
fname = 'AllTagged';
meta = [];

% meta = allSessionMeta(meta,datapth);
meta1 = ALM_SessionMeta(meta,datapth);
meta2 = tjM1_SessionMeta(meta,datapth);
meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)


%% LOAD DATA

for isess = 1:numel(meta)
    disp(' ')
    disp(['Session ' num2str(isess) '/' num2str(numel(meta))])
    disp(' ')
    [sessobj,sesspar(isess)] = loadSessionData(meta(isess),params);
    tag_ = getTagFromObj(sessobj,sesspar(isess),meta(isess));

    % tag(isess) = rmfield(tag_,{'trialid','eventTimes','tagix_allprobes'});
    tag(isess) = tag_;
end


%% SAVE TAG

today = datestr(now, 'yyyymmdd');
fname = [fname '_' params.alignEvent '_' today];
SaveResults(savepath,fname,meta,sesspar,tag)





