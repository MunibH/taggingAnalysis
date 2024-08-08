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
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

params.behav_only = 0;

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


%% LOAD DATA and get nSingle / nMulti

onlyObj = true;
ii = 1;
nSingle = nan(numel(meta),1); % (sessions,1)
nMulti = nan(numel(meta),1);
nTagged = nan(numel(meta),1);
for isess = 1:numel(meta)

    disp(' ')
    disp(['Session ' num2str(isess) '/' num2str(numel(meta))])
    disp(' ')
    [sessobj,sesspar] = loadSessionData(meta(isess),params,onlyObj);

    for iprobe = 1:numel(sessobj.metrics)
        nSingle_(iprobe) = sum(ismember({sessobj.metrics{iprobe}.UnitType},'single'));
        nMulti_(iprobe) = sum(ismember({sessobj.metrics{iprobe}.UnitType},'mua'));
        nTagged_(iprobe) = sum([sessobj.tag(:).probenum]==iprobe);
    end
    nSingle(isess) = sum(nSingle_);
    nMulti(isess) = sum(nMulti_);
    nTagged(isess) = sum(nTagged_);
end

disp(['number of single units: ' num2str(sum(nSingle))] )
disp(['number of multi units: ' num2str(sum(nMulti))] )
disp(['number of tagged units: ' num2str(sum(nTagged))] )

%%

close all

textfs = 12;
cols = getColors;

% Create labels
for isess = 1:numel(meta)
    lab{isess} = [meta(isess).anm ' ' strrep(meta(isess).date,'-','')];
end

% create data
data = [nSingle  nMulti nTagged];

% Create the bar plot with stacked bars
f = figure;
f.Position = [474         365        1108         429];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hBar = bar(ax, data, 'stacked');

% Set the colors
hBar(1).FaceColor = [0.2 0.2 0.2];
hBar(2).FaceColor = [0.6 0.6 0.6];
hBar(3).FaceColor = cols.potent;

hBar(1).EdgeColor = 'none';
hBar(2).EdgeColor = 'none';
hBar(3).EdgeColor = 'none';

xs = hBar(1).XData;

tt1 = text(xs-0.3,sum(data,2)+20,num2str(data(:,2)), "FontSize",textfs,'FontWeight','bold');
tt2 = text(xs-0.3,sum(data,2)+60,num2str(data(:,1)), "FontSize",textfs,'FontWeight','bold');
tt3 = text(xs-0.1,sum(data,2)+100,num2str(data(:,3)), "FontSize",textfs,'FontWeight','bold');
tt4 = text(xs-0.3,sum(data,2)+140,num2str(sum(data,2)), "FontSize",textfs,'FontWeight','bold');
for i = 1:numel(tt1)
    tt1(i).Color = hBar(2).FaceColor;
    tt2(i).Color = hBar(1).FaceColor;
    tt3(i).Color = hBar(3).FaceColor;
    tt4(i).Color = [35, 19, 186]./255;
end

% Set the labels and title
set(gca, 'XTick', xs, 'XTickLabel', lab, 'FontSize', 12, 'XTickLabelRotation', 20); % Align x-ticks with bar positions
xlabel('Session','FontSize',12);
ylabel('Number of units','FontSize',12);
ll = legend('Single units', 'Multi units','PT cells');
ll.Box = 'off';
ll.Location = 'northwest';
ll.FontSize = 12;
ylim([0,800])

