%% Load data
alignmentpth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\alignment\lower_move_thresh_6';
movecorrpth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\results\movement_correlation';

load(fullfile(movecorrpth,'pt_move_corr.mat')) % loads movecorr struct

aligncontents = dir(alignmentpth);
aligncontents = aligncontents(3:end); % remove '.' and '..'
tagalign = [];
for i = 1:numel(aligncontents)
    load(fullfile(aligncontents(i).folder,aligncontents(i).name));
    tagalign = [tagalign tagalignment];
end

%% plot
close all

cols = linspecer(numel(movecorr.feats));

c = getColors;

f = figure;
f.Renderer = 'painters';
t = tiledlayout('flow');
for i = 1:numel(movecorr.feats)
    ax = prettifyAxis(nexttile,2,15);
    hold on;
    feat = movecorr.feats{i};
    f = movecorr.data(i,:);
    temp = corrcoef(f,tagalign);
    cc(i) = temp(1,2);
    scatter(tagalign,f,30,'filled','markeredgecolor','none','markerfacecolor',cols(i,:))
    % title(num2str(round(cc(i),2)))
    axis(ax,'square')
    plot([-1 1],[0 0],'--k')
    plot([0 0],[-1 1],'--k')
    
    p1 = patch([-1 0 0 -1],[-1 -1 1 1],'k');
    p1.EdgeColor = 'none';
    p1.FaceColor = c.potent;
    p1.FaceAlpha = 0.15;
    p2 = patch([0 1 1 0],[-1 -1 1 1],'k');
    p2.EdgeColor = 'none';
    p2.FaceColor = c.null;
    p2.FaceAlpha = 0.15;

    ylabel(strrep(feat,'_','-'),'fontsize',12)

    uistack(p1, 'bottom');
    uistack(p2, 'bottom');


end

xlabel(t,'subspace alignment')









