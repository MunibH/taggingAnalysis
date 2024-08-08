clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - tprime opto out
% - time warp
% - handle tagged units post loading data

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

params.behav_only = 1;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth); % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% subset meta (TODO)

% meta = subsetMetaByParams(meta,params);


%% LOAD DATA

f = figure;
f.Position = [1          41        1920         963];
f.Renderer = 'painters';
t = tiledlayout('flow');

for isess = 1:numel(meta)
    thismeta = meta(isess);
    
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    
    me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
    kin = getKinematics(sessobj, me, sesspar);


    ax = prettifyAxis(nexttile);
    hold on;
    hh = histogram(cell2mat(sessobj.me'));
    hh.EdgeColor = 'none';
    hh.FaceAlpha = 0.7;
    plot([me.moveThresh me.moveThresh],ax.YLim,'Color','r','LineWidth',2 );
    title([thismeta.anm ' ' thismeta.date])
    drawnow;
end
xlabel(t,'Motion energy')
ylabel(t,'Count')










