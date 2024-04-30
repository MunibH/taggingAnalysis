function objvis(meta,params,obj,tag,kin)

% to change rng for trial selection, see ov_balanceAndSplitTrials.m

% TODO
% % clims for kinematics not working...
% % implement time callback

%% Setup GUI
h.bcol = [0.7 0.7 0.7];
h.guicol = [0.7020    0.8549    1.0000];
h.obj = obj;
h.tag = tag;
h.meta = meta;
h.params = params;
h.kin = kin;

clearvars -except h

%% CONTROLS FIGURE
h.fig.control = figure;
set(h.fig.control, 'Units', 'Pixels', 'Position', [15    50   460   941], 'Color', h.bcol);


% conditions to plot
% format: cond  r g b  balanceGroup enabled
tabdat = generateTableData(h.params.condLabel);
h.cond_table = uitable(h.fig.control, 'Data' ,tabdat,'ColumnWidth',{150,40,40,40,60,60}, 'Position', [7,714,450,225], ...
    'ColumnName', {'Cond','Red', 'Green', 'Blue', 'BalanceGroup', 'Enable'}, 'ColumnEditable', true, ...
    'CellEditCallback', {@TableChange, h.fig.control}, 'ColumnFormat', ({[],[],[],[]}),'FontSize',12);


% select units
unitListString = makeUnitListString(h);
h.unit_text = uicontrol(h.fig.control,"Style",'edit','Units', 'Pixels', 'Position', [18,665,163,32], 'String', 'Units', ...
    'FontSize', 10, 'BackgroundColor', h.guicol, 'ForegroundColor','k','FontWeight','bold');
h.unit_list = uicontrol(h.fig.control, 'Style', 'listbox', 'Units', 'Pixels', 'Position', ...
    [18 340 163 325], 'String', unitListString , 'Value', 1, 'BackgroundColor', [1 1 1], ...
    'Max', 3, 'Min', 1, 'Callback',  {@UnitSelect, h.fig.control},'FontSize',9);


% time
tmin = h.params.tmin;
tmax = h.params.tmax;
time_string = ['[' num2str(tmin) ' ' num2str(tmax) ']'];
uicontrol('Style', 'edit', 'Units', 'Pixels', 'Position', [18,290,163,32], 'String', 'Time', ...
    'FontSize', 10, 'BackgroundColor', h.guicol,'ForegroundColor','k','FontWeight','bold');
h.time_edit = uicontrol('Style', 'edit', 'Units', 'Pixels', 'Position', [18,255,163,32], 'String', time_string, ...
    'FontSize', 12, 'BackgroundColor', [1 1 1]);


% smoothing for psth
uicontrol('Style', 'edit', 'Units', 'Pixels', 'Position', [18,210,163,32], 'String', 'Smoothing', ...
    'FontSize', 10, 'BackgroundColor', h.guicol,'ForegroundColor','k','FontWeight','bold');
h.smoothing_edit = uicontrol('Style', 'Edit', 'Units', 'Pixels', 'Position', ...
    [18,175,163,32], 'String', 15, 'Callback', {@updateAxes, h.fig.control},'FontSize',12);


% kinematic feature
uicontrol(h.fig.control,"Style",'edit','Units', 'Pixels', 'Position', [200,665,250,32], 'String', 'KinFeat', ...
    'FontSize', 10, 'BackgroundColor', h.guicol, 'ForegroundColor','k','FontWeight','bold');
h.kin_list = uicontrol(h.fig.control, 'Style', 'listbox', 'Units', 'Pixels', 'Position', ...
    [200 340 250 325], 'String', h.kin.featLeg , 'Value', 1, 'BackgroundColor', [1 1 1], ...
    'Max', 3, 'Min', 1, 'Callback',  {@updateAxes, h.fig.control},'FontSize',12);
uicontrol(h.fig.control,"Style",'edit','Units', 'Pixels', 'Position', [200,290,250,32], 'String', 'KinCLim', ...
    'FontSize', 10, 'BackgroundColor', h.guicol, 'ForegroundColor','k','FontWeight','bold');
h.kin_clim = uicontrol(h.fig.control,"Style",'edit','Units', 'Pixels', 'Position', [200,255,250,32], 'String', '', ...
    'FontSize', 12, 'Callback',  {@updateAxes, h.fig.control});

% save fig button
h.saveFigButton = uicontrol(h.fig.control, 'Style', 'pushbutton', 'Units', 'Pixels', 'Position', [18,50,163,32], 'String', 'SaveFig', ...
    'Callback', {@guisavefig, h.fig.control}, 'FontSize', 12, 'BackgroundColor',h.guicol,'ForegroundColor','k');


%% Init GUI
guidata(h.fig.control, h);

InitPlots(h.fig.control);

% display meta
h = guidata(h.fig.control);
text(h.ax.meta, 0, 0, [h.meta.anm '  ' h.meta.date '  ' strrep(h.meta.region,'_','-')], 'FontSize', 20);

updateAxes([],[],h.fig.control)

end % objvis