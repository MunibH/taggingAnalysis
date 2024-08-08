function timeWarpGUI(~, ~, fig)

h = guidata(fig);
p.parentfig = fig;

p.fig = figure(535);
set(p.fig, 'Units', 'Normalized', 'Position', [0.38 0.5767    0.2194    0.2211], 'Color', h.bcol);

uicontrol('Style', 'Text', 'Units', 'normalized', 'Position', [0.131012658227848,0.625125628140704,0.761392405063291,0.2],...
    'String', '# Licks to Warp (Post Go Cue Licks):', ...
    'FontSize', 14, 'BackgroundColor', h.bcol);

p.nLicks = uicontrol('Style', 'Edit', 'Units', 'normalized', 'Position', ...
    [0.318987341772152,0.465577889447236,0.35,0.2], 'String', 20);


% EXECUTE
uicontrol(p.fig, 'Style', 'pushbutton', 'Units', 'Normalized', 'Position', ...
    [0.282156688959341,0.084924623115578,0.42670407053433,0.206532663316583], 'String', 'Warp Time', ...
    'Callback', {@timeWarp, p.fig(1)}, 'FontSize', 15, 'FontWeight', 'Bold');


guidata(p.fig, p);






end % timeWarpGUI