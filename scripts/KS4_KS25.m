clear,clc,close all

datapth = 'F:\Experiments\';
anm = 'MAH23'; % only MAH23 has ks25 files
dates = {'2024-07-03','2024-07-04','2024-07-05','2024-07-06'};

ii = 1;

for idate = 1:numel(dates)
    thisdate = dates{idate};
    disp(['~~~~~~~~~~ ' thisdate ' ~~~~~~~~~~'])
    [ks4,ks25] = getMetrics(datapth,anm,thisdate);
    for iprobe = 1:numel(ks4)
        disp(['~~~~~~~~~~ ' num2str(ii) ' ~~~~~~~~~~'])
        [rez(ii).ks4,rez(ii).ks25] = getKSMetrics(ks4{iprobe},ks25{iprobe});

        ii = ii + 1;
    end
end

%%
% savepth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\KS4_KS25';
% save(fullfile(savepth,'MAH23_KS4_KS25.mat'),'rez')

%% concatenate results

clear plt

for i = 1:numel(rez)
    plt.ks4.nUnits(i) = rez(i).ks4.nUnits;
    plt.ks4.nSingle(i) = rez(i).ks4.nSingle;
    plt.ks4.nMulti(i) = rez(i).ks4.nMUA;
    plt.ks4.isi{i} = rez(i).ks4.isiv_rate;
    plt.ks4.presenceRatio{i} = rez(i).ks4.presenceRatio;
    plt.ks4.spks_missing_gauss{i} = rez(i).ks4.spks_missing_gauss;

    plt.ks25.nUnits(i) = rez(i).ks25.nUnits;
    plt.ks25.nSingle(i) = rez(i).ks25.nSingle;
    plt.ks25.nMulti(i) = rez(i).ks25.nMUA;
    plt.ks25.isi{i} = rez(i).ks25.isiv_rate;
    plt.ks25.presenceRatio{i} = rez(i).ks25.presenceRatio;
    plt.ks25.spks_missing_gauss{i} = rez(i).ks25.spks_missing_gauss;
end

plt.ks4.isi = cell2mat(plt.ks4.isi');
plt.ks4.presenceRatio = cell2mat(plt.ks4.presenceRatio');
plt.ks4.spks_missing_gauss = cell2mat(plt.ks4.spks_missing_gauss');

plt.ks25.isi = cell2mat(plt.ks25.isi');
plt.ks25.presenceRatio = cell2mat(plt.ks25.presenceRatio');
plt.ks25.spks_missing_gauss = cell2mat(plt.ks25.spks_missing_gauss');

%% plot
close all

c.ks25 = [230, 162, 94]./255;
c.ks4 = [147, 81, 194]./255;
cc(1,:) = c.ks25;
cc(2,:) = c.ks4;

fa = 0.6; % face alpha histogram

% plot nUnits
f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
f.Position = [423         343        1077         591];

% nUnits
ax = prettifyAxis(nexttile);
hold on;
xs = [1 2];
data = {plt.ks25.nUnits,plt.ks4.nUnits};
for i = 1:numel(data)
    this = data{i};
    % b(i) = bar(xs(i),nanmean(this));
    % b(i).FaceColor = clrs(i,:);
    % b(i).EdgeColor = 'none';
    % b(i).FaceAlpha = 1;
    % b(i).BarWidth = 0.7;

    % xx = simple_violin_scatter(xs(i)*ones(size(this)), this, 4, 0.0001);
    xx = xs(i)*ones(size(this));
    scatter(xx, this, 40,'filled', 'markerfacecolor',cc(i,:), 'markeredgecolor',[0.8 0.8 0.8])
    xs_(:,i) = xx;
    ys_(:,i) = this;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end

ax.XTick = xs;
xticklabels({'ks25','ks4'})
ylabel('# of units')
xlim([0.5,2.5])

% nSingle
ax = prettifyAxis(nexttile);
hold on;
xs = [1 2];
data = {plt.ks25.nSingle,plt.ks4.nSingle};
for i = 1:numel(data)
    this = data{i};
    xx = xs(i)*ones(size(this));
    scatter(xx, this, 40,'filled', 'markerfacecolor',cc(i,:), 'markeredgecolor',[0.8 0.8 0.8])
    xs_(:,i) = xx;
    ys_(:,i) = this;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end

ax.XTick = xs;
xticklabels({'ks25','ks4'})
ylabel('# of single units')
xlim([0.5,2.5])

% nMUA
ax = prettifyAxis(nexttile);
hold on;
xs = [1 2];
data = {plt.ks25.nMulti,plt.ks4.nMulti};
for i = 1:numel(data)
    this = data{i};
    xx = xs(i)*ones(size(this));
    scatter(xx, this, 40,'filled', 'markerfacecolor',cc(i,:), 'markeredgecolor',[0.8 0.8 0.8])
    xs_(:,i) = xx;
    ys_(:,i) = this;
end

for i = 1:size(xs_,1)
    patchline(xs_(i,:),ys_(i,:),'EdgeAlpha',0.4)
end

ax.XTick = xs;
xticklabels({'ks25','ks4'})
ylabel('# of mua units')
xlim([0.5,2.5])


% isi
ax = prettifyAxis(nexttile);
hold on;
a = histogram(plt.ks25.isi);
a.FaceAlpha = fa; a.EdgeColor = 'none'; a.FaceColor = c.ks25;
b = histogram(plt.ks4.isi);
b.FaceAlpha = fa; b.EdgeColor = 'none'; b.FaceColor = c.ks4;
xlabel('ISI violation rate (%)')
ylabel('Unit count')
legend({'ks25','ks4'},'location','best','box','off')

% spks missing
ax = prettifyAxis(nexttile);
hold on;
a = histogram(plt.ks25.spks_missing_gauss);
a.FaceAlpha = fa; a.EdgeColor = 'none'; a.FaceColor = c.ks25;
b = histogram(plt.ks4.spks_missing_gauss);
b.FaceAlpha = fa; b.EdgeColor = 'none'; b.FaceColor = c.ks4;
xlabel('Spikes missing (%)')
ylabel('Unit count')

% presence ratio
ax = prettifyAxis(nexttile);
hold on;
a = histogram(plt.ks25.presenceRatio);
a.FaceAlpha = fa; a.EdgeColor = 'none'; a.FaceColor = c.ks25;
b = histogram(plt.ks4.presenceRatio);
b.FaceAlpha = fa; b.EdgeColor = 'none'; b.FaceColor = c.ks4;
xlabel('Presence ratio')
ylabel('Unit count')

fpth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\figs\KS4_KS25';
fn = 'MAH23_KS4_KS25';
mysavefig(f,fpth,fn)

%% Helper functions


function [ks4,ks25] = getMetrics(datapth,anm,date)
ks4fn = ['data_structure_' anm '_' date '.mat'];
if strcmp(anm,'JPV13')
    ks25fn = ['data_structure_' anm '_' date '_ks25.mat'];
else
    ks25fn = ['data_structure_' anm '_' date '_KS25.mat'];
end

% load ks4
load(fullfile(datapth,anm,'Analysis',date,ks4fn))
ks4 = obj.metrics;
clear obj
% load ks25
load(fullfile(datapth,anm,'Analysis',date,ks25fn))
ks25 = obj.metrics;
clear obj

end % getMetrics


function [ks4rez,ks25rez] = getKSMetrics(ks4,ks25)

% num units
ks4rez.nUnits = numel(ks4);
ks25rez.nUnits = numel(ks25);

% isi viol rate
ks4rez.isiv_rate = [ks4.RPV_tauR_estimate]';
ks25rez.isiv_rate = [ks25.RPV_tauR_estimate]';

% % spks missing gaussian
% percentageSpikesMissing_gaussian
ks4rez.spks_missing_gauss = [ks4.percentageSpikesMissing_gaussian]';
ks25rez.spks_missing_gauss = [ks25.percentageSpikesMissing_gaussian]';

% presence ratio
% presenceRatio
ks4rez.presenceRatio = [ks4.presenceRatio]';
ks25rez.presenceRatio = [ks25.presenceRatio]';

% num single and num mua units
% UnitType
ks4rez.nSingle = sum(ismember({ks4.UnitType}','single'));
ks25rez.nSingle = sum(ismember({ks25.UnitType}','single'));
ks4rez.nMUA = sum(ismember({ks4.UnitType}','mua'));
ks25rez.nMUA = sum(ismember({ks25.UnitType}','mua'));


end







