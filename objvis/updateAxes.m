function updateAxes(~, ~, fig)
h = guidata(fig);

% which unit
h.unit = h.unit_list.Value;

% which conditions
[nPossibleConds,~] = size(h.cond_table.Data);
h.cond2plot = find(cell2mat({h.cond_table.Data{:,end}}));

% colors per cond
allcols = cell2mat(reshape({h.cond_table.Data{:,2:4}},nPossibleConds,3));
h.cols = allcols(h.cond2plot,:);

% balance groups
h.balancegroup = cell2mat({h.cond_table.Data{h.cond2plot,5}});

% kinematic feature
h.kinstring = h.kin_list.String{h.kin_list.Value};
h.kinix = h.kin_list.Value;

guidata(fig, h);

updatePSTH([],[],fig);
h = guidata(fig);
updateSpikeRaster([],[],fig);
h = guidata(fig);
updateLickRaster([],[],fig);
h = guidata(fig);
updateKinematics([],[],fig);
h = guidata(fig);
updateISI([],[],fig);
h = guidata(fig);
updateWav([],[],fig);

% % link axes (commenting out for now, making x limits weird)
% linkaxes([h.ax.raster h.ax.psth h.ax.lick h.ax.kin],'x')

% plot event times and set xlims
fields = {'psth','raster','lick','kin'};
xlims = eval(h.time_edit.String);
for i = 1:numel(fields)
    f = fields{i};
    xlim(h.ax.(f),xlims); % xlimits
    plotEventTimes(h.ax.(f),h.params.eventTimes);
end

end







