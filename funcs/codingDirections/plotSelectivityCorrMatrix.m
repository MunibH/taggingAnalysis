function plotSelectivityCorrMatrix(obj,params,dat)


sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay  = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

f = figure; hold on;
imagesc(obj(1).time,obj(1).time,dat);
colorbar; 
% caxis([0 max(max(dat))]);

lw = 1;
ls = '--';
col = [1 1 1] ./ 255;
xline(sample,ls,'Color',col,'LineWidth',lw); yline(sample,ls,'Color',col,'LineWidth',lw)
xline(delay,ls,'Color',col,'LineWidth',lw); yline(delay,ls,'Color',col,'LineWidth',lw)
xline(0,ls,'Color',col,'LineWidth',lw); yline(0,ls,'Color',col,'LineWidth',lw)

xlabel(['Time from ' params.alignEvent ' (s)'])
ylabel(['Time from ' params.alignEvent ' (s)'])
ax = gca;
hold off
colormap(linspecer)
a = colorbar;
a.Label.String = 'Correlation';
ax.FontSize = 10;

axis(ax,'image')



end