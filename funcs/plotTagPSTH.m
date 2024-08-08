% psth
for i = 1:tag.nTag
    [f,ax] = prettyPlot();
    f.Position = [680   560   258   318];
    axRaster = subplot(2,1,1);
    axRaster = prettifyAxis(axRaster);
    axPSTH = subplot(2,1,2);
    axPSTH = prettifyAxis(axPSTH);
    hold(axPSTH,'on');
    hold(axRaster,'on');

    tagdat = mySmooth(squeeze(obj.psth(:,tag.cluid.obj(i),par.cond)),par.sm);

    lasttrial = 0;
    for j = 1:numel(par.cond)
        plot(axPSTH,obj.time,tagdat(:,j),'color',col{j},'linewidth',2)

        if j > 1; pad = par.rasterCondPad; else; pad = 0; end

        % trix = params.trialid{par.cond(j)};
        trix = trials{j};
        thisclunum = tag.cluid.clu(i);
        thisclu = obj.clu{1}(thisclunum);
        mask = ismember(thisclu.trial,trix);
        trial = thisclu.trial(mask);
        trialtm = thisclu.trialtm(mask) - obj.bp.ev.(params.alignEvent)(trial);
        trial = renum(trial) + lasttrial + pad;


        lasttrial = trial(end);

        plot(axRaster,trialtm,trial,'.','color',col{j});

    end

    plotEventTimes(axRaster,params.eventTimes);
    plotEventTimes(axPSTH,params.eventTimes);
    axRaster.FontSize = 11;
    axPSTH.FontSize = 11;
    xlim(axRaster,par.xlims)
    xlim(axPSTH,par.xlims)
    ylabel(axRaster,'Trial')
    ylabel(axPSTH,'Firing rate (spks/s)')
    xlabel(axPSTH,'Time from go cue (s)')
    region = thismeta.region;
    if contains(region,'ALM')
        c = [4, 143, 19]./255;
    elseif contains(region,'tjM1')
        c = [13, 1, 122]./255;
    else
        region = randsample({'L_tjM1','R_tjM1'},1);
        region = region{1};
        c = [13, 1, 122]./255;
    end

    title(axRaster,[thismeta.anm '' thismeta.date '' region ', Unit' num2str(thisclunum)],'fontsize',10,'FontWeight','normal','Color',c,'Interpreter','none');

    if par.sav
        pth = fullfile(utilspth, 'figs', 'TaggedUnits');
        fn = [thismeta.anm '_' thismeta.date '_Unit_' num2str(thisclunum) '_' region '_RasterPSTH'];
        mysavefig(f,pth,fn)
    end


end