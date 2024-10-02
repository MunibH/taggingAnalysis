clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% TODO
% - label licks with which lick number in bout
% - find signifcantly modulated cells during licking/one lick...
%   --- maybe better, find preferred phase for each unit, then sort by that
%       in heatmap

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua'};

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


%% save results (see below for plotting some stuff)
savfig = 1;
savrez = 1;

phasedt = 0.5; % radians

p.vidFs = 400;

p.presampMs = 125; % number ms to extract prior to each lick peak
p.postsampMs = 125; % number ms to extract post to each lick peak
p.presamp = round(p.presampMs/1000 * p.vidFs);
p.postsamp = round(p.presampMs/1000 * p.vidFs);
p.peri_center = 15; % window to look around to correct lick peak

% Parameters
p.filt.Fs = p.vidFs;  % Sampling frequency (in Hz, replace with your actual sampling rate)
p.filt.Fc1 = 0.1;    % Lower cutoff frequency (in Hz)
p.filt.Fc2 = 15;   % Upper cutoff frequency (in Hz)

p.fillmiss_window = 10;

nPerFig = 16;
fpth = 'figs/LickPhaseTuning';
confidence = 0.95;

rpth = 'results/LickPhaseTuning';

%% perform analysis
for isess = 1:numel(meta)
clearvars -except meta utilspth datapth params savfig savrez phasedt p nPerFig fpth confidence isess rpth

    thismeta = meta(isess);

    [obj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(obj,sesspar,thismeta);

    me = loadMotionEnergy(obj, thismeta, sesspar, datapth);
    kin = getKinematics(obj, me, sesspar);


    %% get lick pos, vel, phase

    [lick] = GetLickPhase(obj,p);
    spk = GetUnitSpkPhase(obj,sesspar.cluid,lick);

    %% get binned spiking data around licks and preferred phase    

    for iprobe = 1:numel(spk)
        for iunit = 1:numel(spk{iprobe})
            % (time from lick start, licknum, unit)
            [seq{iprobe}(:,:,iunit),seq_phase] = getSeqLickTriggered(...
                spk{iprobe}(iunit).phase,spk{iprobe}(iunit).licknum,lick.nLicks,phasedt);

            % Assume `phases` is a 1D array containing the phases (in radians) at which a neuron fired spikes
            pref_phase{iprobe}(iunit) = mod(angle(mean(exp(1i * spk{iprobe}(iunit).phase))), 2*pi);

            % test for circular unimodality
            %  [pval z] = circ_rtest(spk{iprobe}(iunit).phase); % for unimodal data
            [otest_pval{iprobe}(iunit)] = circ_otest(sort(spk{iprobe}(iunit).phase)); % for unimodal and axially biomdal date

        end
    end

    %% plots - phase tuning for each unit
    close all

    fname = [thismeta.anm '_' thismeta.date];

    qual = sesspar.quality;

    nClus = sum(cell2mat(cellfun(@numel,sesspar.cluid,'uni',0)));
    nFigs = floor(nClus/nPerFig) + double(rem(nClus,nPerFig)>0);

    f = figure;
    f.Position = [200          84        1161         829];
    f.Renderer = 'painters';
    t = tiledlayout(4,4);
    title(t,strrep(fname,'_',' '))
    figct = 1;
    plotct = 1;
    for iprobe = 1:numel(spk)
        for iunit = 1:numel(spk{iprobe})
            if plotct == nPerFig
                if savfig
                    mysavefig(f,fpth,[fname '_' num2str(figct)]);
                end
                figct = figct + 1;
                f = figure;
                f.Position = [200          84        1161         829];
                f.Renderer = 'painters';
                t = tiledlayout(4,4);
                title(t,strrep(fname,'_',' '))
                plotct = 1;
            end
            ax = nexttile();
            ax.Visible = 'off';
            axPolar = polaraxes(t);
            axPolar.FontSize = 9;
            axPolar.ThetaAxisUnits = 'radians';
            axPolar.Layout.Tile = plotct; % tile location


            phist = polarhistogram(axPolar,spk{iprobe}(iunit).phase,15,'Normalization','probability');
            % Hide the radial axis (r) and grid
            % axPolar.RAxis.Visible = 'off';
            axPolar.ThetaColor = 'k';       % Make the theta (angle) labels black
            % axPolar.RColor = 'none';        % Hide the radial grid lines
            axPolar.ThetaTick = 0:45:360;   % Set the angle ticks (in degrees)
            
            nspks = num2str(numel(spk{iprobe}(iunit).phase));
            title(axPolar,[strrep(thismeta.region{iprobe},'_',' ') ', U' num2str(sesspar.cluid{iprobe}(iunit)) ...
                ', p=' num2str(round(otest_pval{iprobe}(iunit),3)) ', ' qual{iprobe}{iunit} ', ' nspks ' spks'],...
                'fontweight','normal','fontsize',9,'Interpreter','none')

            plotct = plotct + 1;
            drawnow;
        end
    end

    %% save results

    if savrez
        rname = [thismeta.anm '_' thismeta.date];
        SaveResults(rpth, rname, lick,spk,params,sesspar,thismeta)
    end


end






