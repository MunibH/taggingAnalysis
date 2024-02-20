function plotCDProj(rez,obj,params,plt)


sample = mode(obj.bp.ev.sample - obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay - obj.bp.ev.(params.alignEvent));
gc = mode(obj.bp.ev.goCue - obj.bp.ev.(params.alignEvent));


f = figure;
% f.Position = [698   436   343   230];
t = tiledlayout('flow');
for i = 1:numel(rez.cd_labels) % for each coding direction
    ax = nexttile; hold on;
    
    for c = 1:numel(plt.cond2plot)
        trix = params.trialid{plt.cond2plot(c)};
        tempdat = rez.cd_trialdat(:,trix,i); 
        tempdat = mySmooth(tempdat,plt.sm,plt.smtype);
        mu = nanmean(tempdat,2);
        std_err = nanstd(tempdat,[],2) ./ sqrt(size(tempdat,2));
        shadedErrorBar(obj.time,mu,std_err,{'Color',plt.col{c},'LineWidth',plt.lw},plt.alph, ax)
    end

    
    xlim([obj.time(5);2])

    title(rez.cd_labels{i},'FontSize',8)

    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xline(gc,'k--','LineWidth',1)


    curmodename = rez.cd_labels{i};
    if ~strcmpi(curmodename,'ramping')
        shadetimes = obj.time(rez.cd_times.(curmodename));
        x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
        %     y = [-60 -60 50 50];
        fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
        fl.FaceAlpha = 0.3;
        fl.EdgeColor = 'none';
        ylim([y(1) y(3)]);
    end


    hold off

end

xlabel(t,['Time from ' params.alignEvent ' (s)'])
ylabel(t,'Activity (a.u.)')


end