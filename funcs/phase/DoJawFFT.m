function [freq,Y] = DoJawFFT(obj,viddt)

jawpos = [];
for i = 1:obj.bp.Ntrials
    jaw = obj.traj{1}(i).ts(:,2,4);
    jawpos = cat(1,jawpos,jaw);
end

jawpos = fillmissing(jawpos,'previous');

f = figure;
ax = prettifyAxis(gca);
hold on;
Y = fft(jawpos);
freq = (0:length(Y)-1) * (viddt / length(Y));
plot(freq(1:10:end), abs(Y(1:10:end)),'LineWidth',1);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 34.6]);
ylim([0 1.5e5])

end