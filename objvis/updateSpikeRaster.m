function updateSpikeRaster(~,~,fig)

h = guidata(fig);

cla(h.ax.raster)
ylim(h.ax.raster,"auto");

hold(h.ax.raster, 'on');

tm = h.obj.time;

trialid = h.trialid; % these conditions come from updatePSTH
unit = h.unit;
cols = h.cols;

thisclunum = h.params.cluid(unit);
thisclu = h.obj.clu{1}(thisclunum);

rasterCondPad = 5;

lasttrial = 0;
for i = 1:numel(trialid)

    if i > 1; pad = rasterCondPad; else; pad = 0; end

    trix = trialid{i};
    
    mask = ismember(thisclu.trial,trix);
    if sum(mask)==0 % no spikes
        continue
    end
    trial = thisclu.trial(mask);
    trialtm = thisclu.trialtm_aligned(mask);

    trial = renum(trial) + lasttrial + pad;

    lasttrial = trial(end);

    plot(h.ax.raster,trialtm,trial,'.','color',cols(i,:));

end
ylim(h.ax.raster,[h.ax.raster.YLim(1) lasttrial+5])
xlabel(h.ax.raster, ['Time from ' h.params.alignEvent ' (s)']);
ylabel(h.ax.raster, 'Trials');
title(h.ax.raster, 'Spike raster');

end


