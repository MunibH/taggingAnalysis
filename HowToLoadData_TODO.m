clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - time warp
% - fix jaw phase
%   - it now is a cell array, need to interp the time like ME


% for early lick trial formatingg
% - the video frame times for 1 or 2 sessions are a little wonky. see n1
% and n2 below.
%               'JPV11'
%          date: '2023-06-21' 
% above session, side cam for all trials starts at a time that is greater
% than 0.....
% SOLUTION: figure out if those frame times are delayed w.r.t to trial
% start or if they just dropped frames from trial start to then
% then, write a function that identifies this and aligns neural data to
% video data properly. don't hard code that 0.5 anymore that could mean...

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

params.behav_only = 0;

params.useEarly = 1;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);
meta1 = ALM_SessionMeta(meta,datapth);
meta2 = tjM1_SessionMeta(meta,datapth);
meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth);  % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% LOAD DATA

% thismeta = meta(10);

n1 = nan(500,numel(meta));
n2 = nan(500,numel(meta));
for isess = 1:numel(meta)
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params,1);
    sessobj = deleteTaggingTrials(sessobj);
    for i = 1:sessobj.bp.Ntrials
        n1(i,isess) = sessobj.traj{1}(i).frameTimes(1);
        n2(i,isess) = sessobj.traj{2}(1).frameTimes(1);
    end
end

%%

figure; plot(n1-0.5,'.','MarkerSize',10); xlabel('Trial'); ylabel('First frame time (s)'); title('Side cam')
figure; plot(n2-0.5,'.','MarkerSize',10); xlabel('Trial'); ylabel('First frame time (s)'); title('Bottom cam')
[r,c] = find((n1-0.5)>0.05);


%%
tag = getTagFromObj(sessobj,sesspar,thismeta);

me = loadMotionEnergy(sessobj, thismeta, sesspar, datapth);
kin = getKinematics(sessobj, me, sesspar);



%% plot psth

close all

probe = 1;
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
    title(ax1,[sesspar.region{probe}{1} ' U' num2str(sesspar.cluid{probe}(i)) ', ' quality{1}],'Interpreter','none')
   
    plot(ax2,wv,'k')
    plot(ax2,mean(wv,2),'m')

    pause
end
   















