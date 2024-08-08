function CondChanged(~, ~, fig)

h = guidata(fig);
guidata(fig, h);
updateAxes([], [], fig);

end % CondChanged