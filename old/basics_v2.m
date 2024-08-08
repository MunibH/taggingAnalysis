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

params.region = 'any'; % 'alm','tjm1','mc', 'any'
params.probeType = 'any'; % 'h2','np2','np1', 'any'


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
meta = loadMAH23(meta,datapth); % 4 sessions
meta = meta(1);

% subset meta

% meta = subsetMetaByParams(meta,params);

%%  PLOT LICK RASTER and KINEMATICS and Tagged RASTER/PSTH
close all

par.sav = 0; % save figure to .fig and .png
par.xlims = [-2 2];
par.cond = [2,3];
par.sm = 15;
par.rasterCondPad = 10;

par.plotLickRaster = 1;
par.plotPSTH = 0;
par.plotKin = 1;

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

% c1 = [255, 87, 225]./255;
% c2 = [0, 208, 105]./255;
c1 = [1,0,0];
c2 = [0,0,1];
cmap_ = createcolormap(256, c1, c2);

par.feats = {'motion_energy','jaw_ydisp_view1','tongue_angle'};


for isess = 1:numel(meta)
    clearvars -except datapth meta params utilspth isess par col cmap_

    thismeta = meta(isess);

    % LOAD DATA
    [obj,params] = loadSessionData(thismeta,params);
    me = loadMotionEnergy(obj, thismeta, params, datapth);
    kin = getKinematics(obj, me, params);

    % TAGGED UNIT META
    tag.nTag = numel(obj.tag);
    tag.cluid.clu = [obj.tag(:).cluid]; % where tagged units are in obj.clu
    tag.cluid.obj = find(ismember(params.cluid,tag.cluid.clu))'; % where in params.cluid, trialdat, psth

    align = mode(obj.bp.ev.(params.alignEvent));

    trials = params.trialid(par.cond);
    alltrials = cell2mat(trials');

    %%
    % LICK RASTER
    if par.plotLickRaster
        plotLickRaster
        drawnow
    end

    % KINEMATIC heatmaps
    if par.plotKin
        plotKinematics
        drawnow
    end
    

    % % RASTER AND PSTH
    if par.plotPSTH
        plotTagPSTH
        drawnow
    end
    %%


end
