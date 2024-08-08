function plotProj3D(rez,plt)


f = figure;
ax = gca;
hold on;
d = plt.dim2plot3;
for i = 1:numel(plt.cond2plot)
    c = plt.cond2plot(i);
    dat = rez.proj(:,d,c);
    dat = mySmooth(dat,plt.sm,plt.smtype);
    
    plot3(dat(:,1),dat(:,2),dat(:,3), 'Color', plt.col{i}, 'LineWidth', plt.lw)
    for j = 1:numel(rez.eventIX)
        ix = rez.eventIX(j);
        plot3(dat(ix,1),dat(ix,2),dat(ix,3), '.', 'MarkerSize', plt.ms,...
              'Color', plt.col{i}./2)
    end
end
xlabel(['Dim ' num2str(d(1))])
ylabel(['Dim ' num2str(d(2))])
zlabel(['Dim ' num2str(d(3))])

ve = sum(rez.ve(d));
title(['%VE = ' num2str(ve)])

view(ax,[20 30])
grid on

ax.FontSize = 12;


end