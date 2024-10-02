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

%% gather data

all_loss = [];
all_ve = [];
ci_95 = [];
all_neuron_ve = [];
n_neurons = [];
for isess = 1:numel(meta)
    clearvars -except isess meta datapth facemapsvdpth params utilspth all_loss all_ve all_neuron_ve n_neurons ci_95

    load(fullfile(facemapsvdpth, [meta(isess).anm '_' meta(isess).date '_FacemapPredictions.mat']))
    all_loss = cat(2,all_loss,epoch_train_loss);
    all_neuron_ve = cat(2,all_neuron_ve,ve_neuron);
    all_ve = cat(2,all_ve,mean(ve_neuron));
    [~,temp] = mean_CI(ve_neuron);
    ci_95 = cat(2,ci_95,temp);
    n_neurons(isess) = numel(ve_neuron);

end

%% plots
close all

% loss
f = figure;
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
plot(all_loss,'-','LineWidth',1,'Color',[0.25 0.25 0.25])
ylabel('Loss')
xlabel('Epoch')

% ve
f = figure;
f.Position = [680   591   862   287];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
x = n_neurons;
y = all_ve;
e = ci_95;
scatter(x,y,40,'filled','MarkerEdgeColor','none','MarkerFaceColor','k')
errorbar(x,y,e,'LineStyle','none','Color','k','LineWidth',0.5)
ll = line(ax.XLim,[mean(all_ve),mean(all_ve)]);
ll.Color = 'r';
ll.LineStyle = '--';
ll.LineWidth = 2;
ylabel('Variance explained')
xlabel('Number of neurons')


% neuron ve
f = figure;
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
histogram(all_neuron_ve)
mu = mean(all_neuron_ve);
ll = line([mu mu],ax.YLim);
ll.LineStyle = '--';
ll.LineWidth = 2;
ll.Color = 'k';
xlim([-0.2 1])
ylabel('Neuron count')
xlabel('Variance explained')









