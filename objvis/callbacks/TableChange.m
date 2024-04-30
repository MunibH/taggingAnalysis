function TableChange(~, ~, fig)

h = guidata(fig);
guidata(fig, h);
updateAxes([], [], fig);

end