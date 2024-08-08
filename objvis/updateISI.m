function updateISI(~,~,fig)

h = guidata(fig);

cla(h.ax.isi)
ylim(h.ax.isi,"auto");

hold(h.ax.isi, 'on');

unit = h.unit;

thisclunum = h.params.cluid(unit);
thisclu = h.obj.clu{1}(thisclunum);

dt = diff(thisclu.tm);
dt = dt(diff(thisclu.trial)==0);
if isempty(dt)
    dt = 1;
end

ISIedges = -0.02:0.0005:0.02;
Nisi = histc([dt; -dt], ISIedges);
ISIcoord = ISIedges+mean(diff(ISIedges))./2;

mx = max(Nisi);
mx = 1.05.*mx;

hold(h.ax.isi, 'off');
thr = 0.0025;
b = bar(h.ax.isi, ISIcoord, Nisi); 
h.ax.isi = prettifyAxis(h.ax.isi,2); % bar undoes prettifyAxis()
hold(h.ax.isi, 'on');
ylims = h.ax.isi.YLim;
f = fill(h.ax.isi,[-thr thr thr -thr], [0 0 ylims(2) ylims(2)], [1 0 0]);
% ylim(ylims)
set(b, 'LineStyle', 'none', 'BarWidth', 1, 'FaceColor', 'b');
set(f, 'Linestyle', 'none');
set(f,'FaceAlpha',0.3)
xlabel(h.ax.isi,'ISI')
ylabel(h.ax.isi,'# spks')


end


