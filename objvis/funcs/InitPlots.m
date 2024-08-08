function InitPlots(fig)
h = guidata(fig);

h.fig.main = figure;

set(h.fig.main, 'Units', 'Normalized', 'Position', [0.2469    0.0565    0.7500    0.8602]);

h.ax.raster = axes;
set(h.ax.raster, 'Units', 'Normalized', 'Position', [0.05 0.095 0.3 0.43], 'Color', [1 1 1]);
h.ax.raster = prettifyAxis(h.ax.raster,2);

h.ax.psth = axes;
set(h.ax.psth, 'Units', 'Normalized', 'Position', [0.05 0.65 0.3 0.3], 'Color', [1 1 1]);
h.ax.psth = prettifyAxis(h.ax.psth,2);

h.ax.lick = axes;
set(h.ax.lick, 'Units', 'Normalized', 'Position', [0.4 0.095 0.3 0.43], 'Color', [1 1 1]);
h.ax.lick = prettifyAxis(h.ax.lick,2);

h.ax.kin = axes;
set(h.ax.kin, 'Units', 'Normalized', 'Position', [0.4 0.65 0.3 0.3], 'Color', [1 1 1]);
h.ax.kin = prettifyAxis(h.ax.kin,2);

h.ax.isi = axes;
set(h.ax.isi, 'Units', 'Normalized', 'Position', [0.75 0.65 0.2 0.28], 'Color', [1 1 1]);
h.ax.isi = prettifyAxis(h.ax.isi,2);

h.ax.wav = axes;
set(h.ax.wav, 'Units', 'Normalized', 'Position', [0.75 0.2 0.2 0.32], 'Color', [1 1 1]);
h.ax.wav = prettifyAxis(h.ax.wav,2);

h.ax.meta = axes;
set(h.ax.meta, 'Units', 'Normalized', 'Position', [0.74 0.09 0.2 0.05], 'Color', [1 1 1], 'Visible','off');

guidata(fig, h);

end