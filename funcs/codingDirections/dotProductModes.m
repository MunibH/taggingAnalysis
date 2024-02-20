function dotProductModes(rez,modes)

modenames = rez.cd_labels;

% dot product b/w modes (measure of orthogonality)
dots = modes' * modes;
f = figure;
imagesc(dots);
colorbar
cmap = flip(gray);
colormap(cmap)
caxis([min(min(dots)),max(max(dots))])
xticks([1:1+numel(modenames)])
xticklabels(modenames)
yticks([1:1+numel(modenames)])
yticklabels(modenames)
title('Orthogonality of CDs')
end % dotProductModes