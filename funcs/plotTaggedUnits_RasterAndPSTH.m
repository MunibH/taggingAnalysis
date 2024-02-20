function plotTaggedUnits_RasterAndPSTH()

%% DATA and PARAMETERS

params = defaultParams();
% % specify changes here
% params.alignEvent = 'firstLick';

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
% meta = loadJPV13(meta,datapth);

% meta = meta(1);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

disp("LOADING ALL DATA OBJECTS")
obj = loadObjs(meta);

%% PROCESS TAGGED UNITS ONLY

prbnum = 1; % I only have single probe recordings in my data,
% you'll have to rewrite some of this code to handle dual probe recordings

% tag will be like obj, but for one entry for each session
% and only contains tagged unit data
% this makes it faster to process all the data if you just want to look at
% tagged units
tag = getTagStruct(meta,obj,params);


%% PLOT
close all

nTagged = numel([tag(:).cluid]);

% event times
eventTimes = getEventTimes(obj(1).bp.ev,params(1).events,params.alignEvent);
align = mode(obj(1).bp.ev.(params.alignEvent));

cond2plot = [2,3]; % rhit,lhit

cols(1,:) = [0,0,1];
cols(2,:) = [1,0,0];

rasterCondPad = 10;
sm = 11;

for isess = 1:numel(tag)
    for iunit = 1:numel(tag(isess).cluid)
        f = figure;
        f.Renderer = 'painters';
        % psth
        axRaster = subplot(2,1,1);
        axPSTH = subplot(2,1,2);
        hold(axPSTH,'on');
        hold(axRaster,'on');

        tagdat = mySmooth(squeeze(tag(isess).psth(:,iunit,cond2plot)),sm);

        lasttrial = 0;
        for j = 1:numel(cond2plot)
            plot(axPSTH,tag(isess).time,tagdat(:,j),'color',cols(j,:),'linewidth',2)

            if j > 1; pad = rasterCondPad; else; pad = 0; end

            trix = tag(isess).trialid{cond2plot(j)};
            thisclu = tag(isess).clu(iunit);
            mask = ismember(thisclu.trial,trix);
            trial = thisclu.trial(mask);
            trialtm = thisclu.trialtm_aligned(mask);% - obj(isess).bp.ev.(params(isess).alignEvent)(trial);
            trial = renum(trial) + lasttrial + pad;


            lasttrial = trial(end);

            plot(axRaster,trialtm,trial,'.','color',cols(j,:));

        end
        plotEventTimes(axRaster,eventTimes);
        plotEventTimes(axPSTH,eventTimes);
        axRaster.FontSize = 11;
        axPSTH.FontSize = 11;
        xlim(axRaster,[-2,2])
        xlim(axPSTH,[-2,2])
        % title(axRaster,['Unit ' num2str(thisclunum)],'fontsize',10)
        ylabel(axRaster,'Trial')
        ylabel(axPSTH,'Firing rate (spks/s)')
        xlabel(axPSTH,'Time from go cue (s)')


    end



end

