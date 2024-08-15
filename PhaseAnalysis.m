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


%% LOAD DATA

thismeta = meta(1);

[obj,sesspar] = loadSessionData(thismeta,params);
tag = getTagFromObj(obj,sesspar,thismeta);

me = loadMotionEnergy(obj, thismeta, sesspar, datapth);
kin = getKinematics(obj, me, sesspar);


%% get licks
close all
clear all_licks lick_trialtm lick_trial

% for each trial
% find lick contact times
% get jaw pos -100 to +100 ms around each lick contact
% calculate phase?

vidFs = 400;
presamp = 50;
postsamp = 50;
presampms = presamp/vidFs * 1000;
postsampms = postsamp/vidFs * 1000;
peri_center = 15;

% Parameters
Fs = vidFs;  % Sampling frequency (in Hz, replace with your actual sampling rate)
Fc1 = 0.1;    % Lower cutoff frequency (in Hz)
Fc2 = 15;   % Upper cutoff frequency (in Hz)

% Normalized cutoff frequencies (between 0 and 1)
Wn = [Fc1 Fc2] / (Fs / 2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(4, Wn, 'bandpass');

traj = obj.traj{1}; % side cam
feat = 'jaw';
ifeat = ismember(traj(1).featNames,feat);

lickct = 1;
for itrial = 1:obj.bp.Ntrials
    licktimes = sort([obj.bp.ev.lickL{itrial} obj.bp.ev.lickR{itrial}]);

    jaw = traj(itrial).ts(:,2,ifeat); % jaw, side cam, ycoord
    ft = traj(itrial).frameTimes - 0.5;

    licktimes(licktimes<ft(1) | licktimes>ft(end)) = [];

    % f = figure;
    % ax = gca;
    % hold on;
    % plot(ft,jaw)
    % for i = 1:numel(licktimes)
    %     plot([licktimes(i) licktimes(i)], ax.YLim, 'k-' )
    % end


    % f = figure;
    % ax = gca;
    % hold on;
    for ilick = 1:numel(licktimes)
        % yyaxis left
        % cla(ax)
        % yyaxis right
        % cla(ax)

        tix = findTimeIX(ft,licktimes(ilick));
        startix = max(1,tix-presamp);
        endix = min(tix+postsamp,numel(ft));
        lick = jaw(startix:endix);

        % find start and end time of 'lick'
        licktimes_ft = ft(startix:endix);


        % bandpass filter (fill missing data first)
        x_filled = fillmissing(lick, 'movmean', 10);  % Replace window_size with an appropriate value
        if any(isnan(x_filled)) % there was one lick I found that fillmissing still didn't fill all values, skip that one...
            continue
        end
        lick_filt = filtfilt(b,a,x_filled);

        % where bpod detected lick
        center_idx = round(numel(lick_filt)/2);
        % correct the center idx (bpod lick time might not actually be peak of lick
        [~,center_idx_] = max(lick_filt(center_idx-peri_center:center_idx+peri_center));
        center_idx_corrected = center_idx_ + center_idx - peri_center;
        % walk bckward from center find start of lick
        [~, min_back_idx] = min(lick_filt(1:center_idx_corrected));  % Index of minimum before the peak
        % walk forward from center find end of lick
        [~, min_forward_idx] = min(lick_filt(center_idx_corrected:end));  % Index of minimum after the peak
        min_forward_idx = min_forward_idx + center_idx_corrected - 1;  % Adjust index relative to the original signal


        % yyaxis left
        % plot(lick)
        % xline(center_idx_corrected)
        % xline(min_back_idx)
        % xline(min_forward_idx)
        % yyaxis right
        % plot(lick_filt)

        all_licks{lickct} = lick_filt(min_back_idx:min_forward_idx);
        ix.lick_center{lickct} = center_idx_corrected;
        ix.lick_center{lickct} = center_idx_corrected;
        ix.lick_center{lickct} = center_idx_corrected;


        % now I have the lick, I need the lick start and end times within
        % the trial (frameTimes -- ft)
        lick_trialtm{lickct} = licktimes_ft(min_back_idx:min_forward_idx);
        lick_trial(lickct) = itrial;
        lickct = lickct + 1;

    end


end

% % resample all licks
% resampledLickData = cell2mat(resampleToMedianLength(all_licks))';



%% plot licks

close all
% cm = cmap_sweep(10,jet);
% cm = jet(10);
cm = flipud(thermal);
cm = cm(1:end,:);
nCols = size(cm,1);
nLicks = numel(all_licks);

% Create a query points vector for interpolation
query_points = linspace(1, nCols, nLicks);

% Interpolate each of the R, G, B channels separately
cm = interp1(1:nCols, cm, query_points);

% Ensure the result is a valid colormap (in case of rounding errors)
% cm = max(0, min(1, cm));  % Ensure values are between 0 and 1

f = figure;
ax = prettifyAxis(gca);
hold on;
for i = 1:nLicks
    tvec = (1:numel(all_licks{i}))./vidFs * 1000;
    patchline(tvec,zscore(all_licks{i}),'EdgeColor',cm(i,:))
end
xlabel('Time (ms)')
ylabel('zscored licks')
xlim([0 250])
% 
% f = figure;
% ax = prettifyAxis(gca);
% hold on;
% % phase = linspace(0,2*pi,size(resampledLickData,1));
% for i = 1:nLicks
%     patchline(phase,zscore(resampledLickData(:,i)),'EdgeColor',cm(i,:))
% end
% xlabel('Phase (rad)')
% ylabel('zscored licks')
% xlim([0 2*pi])

%% plot licks with phase on x-axis

close all

% phase of licks
phases = assignPhases(all_licks);

f = figure;
ax = prettifyAxis(gca);
hold on;
for i = 1:nLicks
    patchline(phases{i},zscore(all_licks{i}),'EdgeColor',cm(i,:))
end
xlabel('Phase (rad)')
ylabel('zscored licks')
xlim([0 2*pi])

%% trigger PT units on start of lick
clear useSpikes nSpikesOnLick spk_licktm spk_licktm_ts lick_ts phase_licktm_ts spk_lick_num spk_lick_tm spk_lick_phase

lick_durations = cell2mat( cellfun(@numel,lick_trialtm,'uni',0)  );
longest_lick = max(lick_durations);
longest_lick_duration = longest_lick*1/vidFs;
lick_time_vec = (0:(longest_lick-1)).*(1/vidFs);

nLicks = numel(lick_trialtm);

nProbes = numel(thismeta.probe);

spk_licktm_ts = cell(nProbes,1);
spk_lick_num = cell(nProbes,1);
spk_lick_tm = cell(nProbes,1);
spk_lick_phase = cell(nProbes,1);
unitct = 1;
for iprobe = 1:nProbes
    nTag = tag.nTag(iprobe);
    spk_licktm_ts{iprobe} = zeros(longest_lick,nLicks,nTag); % (time,lick,unit)
    for itag = 1:nTag
        spk_lick_num{iprobe}{itag} = [];
        spk_lick_tm{iprobe}{itag} = [];
        spk_lick_phase{iprobe}{itag} = [];

        cluid = tag.id.clu{iprobe}(itag);
        clu = obj.clu{iprobe}(cluid);
        for ilick = 1:nLicks
            thistrial = lick_trial(ilick);
            thislick_trialtm = lick_trialtm{ilick};
            % get spikes within lick times on this trial

            onTrial = clu.trial==thistrial;
            inTime = (clu.trialtm>=thislick_trialtm(1)) & (clu.trialtm<=thislick_trialtm(end));

            useSpikes = (onTrial & inTime);
            nSpikesOnLick(ilick,unitct) = sum(useSpikes);

            spk_trialtm = clu.trialtm(useSpikes);

            % align to lick onset
            spk_licktm{iprobe}{itag}{ilick} = spk_trialtm - thislick_trialtm(1);
            spk_lick_num{iprobe}{itag} = cat(1, spk_lick_num{iprobe}{itag}, ones(size(spk_trialtm))*ilick );
            spk_lick_tm{iprobe}{itag} = cat(1, spk_lick_tm{iprobe}{itag}, spk_trialtm );

            for i = 1:numel(spk_licktm{iprobe}{itag}{ilick})
                lix = findTimeIX(lick_time_vec,spk_licktm{iprobe}{itag}{ilick}(i));
                spk_licktm_ts{iprobe}(lix,ilick,itag) = 1;
                spk_lick_phase{iprobe}{itag} = cat(1, spk_lick_phase{iprobe}{itag},phases{ilick}(lix));
            end
        end
        unitct = unitct + 1;
    end
end

% put all licks into a time series
lick_ts = nan(longest_lick,nLicks);
for i = 1:nLicks
    thislick = all_licks{i};
    lick_ts(1:numel(thislick),i) = thislick;
end


%% plots - num spikes per lick

%
f = figure;
f.Position = [879   431   264   251];
f.Renderer = 'painters';
ax = prettifyAxis(gca);
hold on;
xs = 1:size(nSpikesOnLick,2);
for i = 1:numel(xs)
    this = nSpikesOnLick(:,i);
    xx = simple_violin_scatter(xs(i)*ones(size(this)),this,1000,0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor','k', 'markeredgecolor','none')
end
ax.XTick = xs;
xticklabels({'Tag1','Tag2','Tag3','Tag4','Tag5'})
ylabel('# of spikes per lick')
ylim([-1,ax.YLim(2)])

%% plots - lick triggered raster

close all

for iprobe = 1:nProbes
    nTag = tag.nTag(iprobe);
    for itag = 1:nTag
        f = figure;
        f.Position = [1238         438         264         383];
        f.Renderer = 'painters';
        ax = prettifyAxis(gca);
        hold on;
        scatter(spk_lick_phase{iprobe}{itag},spk_lick_num{iprobe}{itag},'k.')

        xlabel('Phase (rad)')
        ylabel('Lick #')
        xlim([0 2*pi])
    end
end

%% plots - mean spike rate lick triggered, about phase
close all
phasedt = 0.05;
sm = 21;

for iprobe = 1:nProbes
    nTag = tag.nTag(iprobe);
    for itag = 1:nTag
        f = figure;
        f.Position = [680   619   368   259];
        f.Renderer = 'painters';
        ax = prettifyAxis(gca);
        hold on;

        [seq,seq_phase] = getSeqLickTriggered(...
            spk_lick_phase{iprobe}{itag},spk_lick_num{iprobe}{itag},phasedt,sm);

        [m,h] = mean_CI(seq);
        shadedErrorBar(seq_phase,m,h,{'Color','k','LineWidth',2},0.3,ax);
        plot([pi pi],ax.YLim,'k--','LineWidth',2)
        xlabel('Phase (rad)')
        ylabel('spks/phase')
        xlim([0 2*pi])
    end
end




%% plots - lick triggered heatmap

% close all
% %
% 
% lick_durations = cell2mat( cellfun(@numel,lick_trialtm,'uni',0)  );
% 
% cm = cividis;
% 
% for iprobe = 1:nProbes
%     nTag = tag.nTag(iprobe);
%     for itag = 1:nTag
%         f = figure;
%         f.Position = [680   379   419   499];
%         f.Renderer = 'painters';
%         ax = prettifyAxis(gca);
%         hold on;
% 
% 
%         this = spk_licktm_ts{iprobe}(:,:,itag);
%         imagesc(lick_time_vec*1000,1:size(this,2),this');
%         colormap(cm)
% 
%         % for ilick = 1:nLicks
%         %     thislick = all_licks{ilick};
%         %     [~,maxix] = max(thislick);
%         %     maxtime = maxix*(1/vidFs);
%         %     plot(maxtime*1000,ilick,'w.','MarkerSize',5)
%         %     lickendtime = numel(thislick)*(1/vidFs);
%         %     plot(lickendtime*1000,ilick,'m.','MarkerSize',3)
%         % end
% 
%         xlabel('Time from lick (ms)')
%         ylabel('Lick #')
%         xlim([0 longest_lick_duration*1000])
%         ylim([-0.1 nLicks+0.1])
%     end
% 
% end


