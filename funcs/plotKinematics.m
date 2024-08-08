[f,ax] = prettyPlot();
f.Position = [398         370        1008         378];
t = tiledlayout('flow');

for ifeat = 1:numel(par.feats)
    feat = par.feats{ifeat};
    featix = find(ismember(kin.featLeg,feat));
    kindat = kin.dat(:,alltrials,featix);
    kindat = removeOutliers(kindat,3);

    ax = nexttile;
    ax = prettifyAxis(ax);
    hold on;


    imagesc(ax,obj(1).time,1:size(kindat,2),kindat','interpolation','bilinear')
    colormap(ax,cmap_)
    c = colorbar; c.Label.String = feat;
    plotEventTimes(ax,params.eventTimes)
    ax.FontSize = 11;
    xlim(par.xlims)
    yl = cumsum(cell2mat(cellfun(@(x) numel(x),trials,'uni',0)));
    for i = 1:numel(trials)-1
        yline(yl(i),'k--')
    end
    ax = prettifyAxis(ax);
    ylim([1,numel(alltrials)])

end
ylabel(t,'Trials')
xlabel(t,['Time from ' params.alignEvent '(s)'])

if par.sav
    pth = fullfile(utilspth, 'figs', thismeta.anm, thismeta.date);
    fn = ['kinematics'];
    mysavefig(f,pth,fn)
end

% KINEMATIC trial-avgs
[f,ax] = prettyPlot();
f.Position = [321   223   922   194];
t = tiledlayout('flow');

for ifeat = 1:numel(par.feats)
    feat = par.feats{ifeat};
    featix = find(ismember(kin(1).featLeg,feat));
    ax = nexttile;
    ax = prettifyAxis(ax);
    hold on;

    for icond = 1:numel(par.cond)
        kindat = kin.dat(:,trials{icond},featix);
        % kindat = removeOutliers(kindat,3);
        kindat = squeeze(nanmean(kindat,2));

        plot(ax,obj(1).time,kindat,'Color',col{icond},'LineWidth',1.5)
    end

    plotEventTimes(ax,params.eventTimes)
    ylabel(par.feats{ifeat},'Interpreter','none')
    ax.FontSize = 11;
    xlim(par.xlims)
end
xlabel(t,['Time from ' params.alignEvent '(s)'])

if par.sav
    pth = fullfile(utilspth, 'figs', thismeta.anm, thismeta.date);
    fn = 'kinematics-tavg';
    mysavefig(f,pth,fn)
end
