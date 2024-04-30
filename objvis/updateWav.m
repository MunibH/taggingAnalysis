function updateWav(~,~,fig)

h = guidata(fig);

cla(h.ax.wav)
ylim(h.ax.wav,"auto");

hold(h.ax.wav, 'on');

unit = h.unit;

thisclunum = h.params.cluid(unit);
thisclu = h.obj.clu{1}(thisclunum);

if ~isfield(thisclu,'spkWavs')
    x = h.ax.wav.XLim(1);
    y = h.ax.wav.YLim(2)/3;
    text(h.ax.wav,x,y,'spike waveforms not found',...
        'FontSize', 20)
    return
end

wavs = squeeze(thisclu.spkWavs);
meanWav = mean(wavs,2);

plot(h.ax.wav,wavs,'Color',[0,0,0,0.3],'LineWidth',1.5)
hold(h.ax.wav,'on')
plot(h.ax.wav,meanWav,'m','LineWidth',3)
hold(h.ax.wav,'off')
title(h.ax.wav,'Waveforms')

end
