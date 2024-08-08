function updateLickRaster(~,~,fig)

h = guidata(fig);

cla(h.ax.lick)
ylim(h.ax.lick,"auto");

hold(h.ax.lick, 'on');

tm = h.obj.time;

trialid = h.trialid; % these conditions come from updatePSTH
cols = h.cols;


rasterCondPad = 5;

align = mode(h.obj.bp.ev.(h.params.alignEvent));

lick.l = cellfun(@(x) x-align, h.obj.bp.ev.lickL,'uni',0);
lick.r = cellfun(@(x) x-align, h.obj.bp.ev.lickR,'uni',0);


lasttrial = 0;
for j = 1:numel(trialid)
    clear allLicks
    if j > 1; pad = rasterCondPad; else; pad = 0; end

    trix = trialid{j};
    lick.trialtm.l = lick.l(trix);
    lick.trialtm.r = lick.r(trix);
    lick.trial = renum(trix) + lasttrial + pad;

    lasttrial = lick.trial(end);

    % need to make lick.trialtm.l/r and lick.trial same size

    allLicks.ltm = cell2mat(lick.trialtm.l');
    allLicks.rtm = cell2mat(lick.trialtm.r');

    rindex = 1;
    lindex = 1;
    for itrial = 1:numel(trix)
        nLicksL = numel(lick.trialtm.l{itrial});
        nLicksR = numel(lick.trialtm.r{itrial});
        allLicks.lTrial(lindex:lindex + nLicksL) = lick.trial(itrial);
        allLicks.rTrial(rindex:rindex + nLicksR) = lick.trial(itrial);
        rindex = rindex + nLicksR;
        lindex = lindex + nLicksL;
    end
    allLicks.rTrial = allLicks.rTrial(1:end-1);
    allLicks.lTrial = allLicks.lTrial(1:end-1);


    plot(h.ax.lick,allLicks.rtm,allLicks.rTrial,'.','color','b');
    plot(h.ax.lick,allLicks.ltm,allLicks.lTrial,'.','color','r');

end

ylim(h.ax.lick,[h.ax.lick.YLim(1) lasttrial+5])
xlabel(h.ax.lick, ['Time from ' h.params.alignEvent ' (s)']);
ylabel(h.ax.lick, 'Trials');
title(h.ax.lick, 'Lick raster');

end


