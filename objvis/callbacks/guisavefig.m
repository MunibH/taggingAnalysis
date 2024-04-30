function guisavefig(~,~,fig)

h = guidata(fig);
f = h.fig.main;

defpth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis\figs\objvis';
nowstring = datestr(now, 'yyyymmdd_HHMMSS');
savePath = fullfile(defpth,[h.meta.anm '_' h.meta.date '_saved_' nowstring]);

[file,path] = uiputfile('*.*','Specify directory and file name',savePath);

saveas(f,fullfile(path,file),'png')


end








