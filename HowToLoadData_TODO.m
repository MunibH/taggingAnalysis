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
% - fix jaw phase
%   - it now is a cell array, need to interp the time like ME

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

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

thismeta = meta(1);

[sessobj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(sessobj,sesspar,thismeta);

me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
kin = getKinematics(sessobj, me, sesspar);



%% plot psth

close all

probe = 2;
nUnits = numel(sesspar.cluid{probe});

cols = getColors;
cond2plot = [2,3];
c(1,:) = cols.rhit;
c(2,:) = cols.lhit;

xl = [-2.1,params.tmax];

f = figure;
f.Position = [370         577        1073         301];
ax1 = prettifyAxis(subplot(1,2,1));
hold on;
ax2 = prettifyAxis(subplot(1,2,2));
hold on;
for i = 1:nUnits
    cla(ax1)
    cla(ax2)
    quality = sesspar.quality{probe}(i);
    wv = squeeze(sessobj.clu{probe}(sesspar.cluid{probe}(i)).spkWavs);
    for j = 1:numel(cond2plot)
        cond = cond2plot(j);
        thispsth = squeeze(sessobj.psth{probe}(:,i,j));
        plot(ax1,sessobj.time,thispsth,'Color',c(j,:),'LineWidth',2)
    end
    plotEventTimes(ax1,sesspar.eventTimes,'k',false)
    xlim(ax1,xl)
    title(ax1,[sesspar.region{probe} ' U' num2str(sesspar.cluid{probe}(i)) ', ' quality{1}],'Interpreter','none')
   
    plot(ax2,wv,'k')
    plot(ax2,mean(wv,2),'m')

    pause
end
   















