function plotProj2D(rez,plt,params)


f = figure;
t = tiledlayout('flow');
for i = 1:numel(plt.dim2plot2)
    d = plt.dim2plot2(i);
    ax = nexttile;
    hold on;
    for c = 1:numel(plt.cond2plot)
        trix = params.trialid{plt.cond2plot(c)};
        temp = squeeze(rez.proj_trialdat(:,trix,i)); % single trial data for all trials in a condition, one cluster
        temp = mySmooth(temp,plt.sm,plt.smtype); % smooth with causal gaussian filter (smtype is optional argument)
        mu = mean(temp,2); % mean over trials in condition
        std_err = std(temp,[],2) ./ sqrt(size(temp,2)); % standard error
        shadedErrorBar(plt.time,mu,std_err,{'LineWidth',plt.lw,'Color',plt.col{c}},plt.alph,ax)
    end
    for j = 1:numel(rez.eventIX)
        xline(plt.time(rez.eventIX(j)),'k--');
    end
    title(['Dim ' num2str(d) ' | %VE=' num2str(rez.ve(d)) ])
end
xlabel(t,['Time from ' params.alignEvent ' (s)'])
ylabel(t, 'Activity (a.u.)')



end