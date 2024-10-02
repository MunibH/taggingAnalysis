%%

addpath(genpath('C:\npy-matlab\npy-matlab'))

datapth = 'F:\Experiments\JPV13\Analysis\2023-10-03\Probe1\kilosort4';
ampfn = 'amplitudes.npy';

amp = readNPY(fullfile(datapth,ampfn));

%%

datapth = 'F:\Experiments\JPV13\Analysis\2023-10-03';
fn = 'Probe1Data.mat';
load(fullfile(datapth,fn))

%% 

clus = unique(ks.spike_clusters);

f = figure;
f.Renderer = 'painters';
f.Position = [486         430        1149         340];
ax = prettifyAxis(gca);
hold on;
for i = 1:numel(clus)
    cla(ax)

    thisclu = clus(i);
    ix = ks.spike_clusters==thisclu;
    times = ks.spike_times(ix);
    thisamp = amp(ix);
    scaling = ks.tempScalingAmps(ix);
    

    plot(times,thisamp.*scaling/3,'.','MarkerSize',10)

    % duration = times(end)-times(1);
    xs = linspace(times(1),times(end),20);
    for j = 1:numel(xs)
        plot([xs(j),xs(j)],ax.YLim,'k')
    end
    
    xlabel('Time (s)')
    ylabel('Spike amplitude (\muV)')
    title(['JPV13, 20231003, Unit ' num2str(i)],'fontweight','normal')
    pause

end