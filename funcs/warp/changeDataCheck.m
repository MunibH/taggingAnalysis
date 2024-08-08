function changeDataCheck(~, ~, fig)
h = guidata(fig);

probe = get(h.probeList, 'Value');
unit = get(h.unitList, 'Value');

clu = h.obj.clu{probe}(unit);

dataToUse = h.psthDataList.String{h.psthDataList.Value};
if ~isfield(h.obj.clu{probe}(unit),dataToUse)
    f = msgbox('Make sure to align or warp data first');
    return
end

updateRaster(fig, clu);
updatePSTH(fig, clu);
updateVideo([],[],fig);
updateISI(fig, clu);
updateWav(fig, clu);
updateBehav([],[],fig);

figure(h.fig(1))






