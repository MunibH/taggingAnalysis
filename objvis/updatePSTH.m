function updatePSTH(~,~,fig)

h = guidata(fig);

cla(h.ax.psth)
ylim(h.ax.psth,"auto");

hold(h.ax.psth, 'on');

tm = h.obj.time;

sm = str2double(h.smoothing_edit.String);

% avg and stderr over equal number of trials across conditions
nConds = numel(h.cond2plot);
trials = ov_balanceAndSplitTrials(h.params.trialid, h.cond2plot, h.balancegroup);
for i = 1:nConds
    this = squeeze(h.obj.trialdat(:,h.unit,trials{i}));
    this = mySmooth(this,sm,'reflect'); 

    mu = mean(this,2);
    stderr = std(this,[],2)./sqrt(size(this,2));
    shadedErrorBar(tm, mu, stderr, {'Color',h.cols(i,:),'LineWidth',2}, 0.3, h.ax.psth);
end

xlabel(h.ax.psth, ['Time from ' h.params.alignEvent ' (s)']);
ylabel(h.ax.psth, 'spks/sec');
title(h.ax.psth, ['Unit ' num2str(h.unit)]);

% plotEventTimes(h.ax.psth,h.params.eventTimes)

h.trialid = trials;

guidata(fig, h);


end



