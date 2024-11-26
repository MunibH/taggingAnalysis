f = figure;
f.Position = [943   363   292   335];
f.Renderer = "painters";
hold on;
ax = prettifyAxis(gca);
hold on;
xs = [1 2 4 5];
for i = 1:numel(xs)
    this = toplot{i};
    xx = simple_violin_scatter(xs(i)*ones(size(this)), this, numel(this)./i, 0.5);
    scatter(xx, this, 8,'filled', 'markerfacecolor',c(i,:), 'markeredgecolor','none')
end
ax.XTick = xs;
xticklabels({'null','nullpt','potent','potentpt'})
ylabel('Subspace contribution')
