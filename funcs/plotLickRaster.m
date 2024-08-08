[f,ax] = prettyPlot();
f.Position = [680   560   258   318];
ax = prettifyAxis(ax);
hold on;

lick.l = cellfun(@(x) x-align, obj.bp.ev.lickL,'uni',0);
lick.r = cellfun(@(x) x-align, obj.bp.ev.lickR,'uni',0);

lasttrial = 0;
for j = 1:numel(par.cond)
    clear allLicks
    if j > 1; pad = par.rasterCondPad; else; pad = 0; end

    % trix = params.trialid{par.cond(j)};
    trix = trials{j};
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


    plot(ax,allLicks.rtm,allLicks.rTrial,'.','color',col{1});
    plot(ax,allLicks.ltm,allLicks.lTrial,'.','color',col{2});

end
plotEventTimes(ax,params.eventTimes);
ax.FontSize = 11;
xlim(ax,par.xlims)
title([thismeta.anm ' ' thismeta.date],'fontsize',10)
ylabel(ax,'Trial')
xlabel(ax,'Time from go cue (s)')

if par.sav
    pth = fullfile(utilspth, 'figs', thismeta.anm, thismeta.date);
    fn = ['lickRaster'];
    mysavefig(f,pth,fn)
end