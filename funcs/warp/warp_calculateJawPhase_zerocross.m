function jawphase = warp_calculateJawPhase_zerocross(dlc,jawix,coord,gocueTimes,lickL)

%%

close all
f = figure;
ax1 = subplot(2,1,1);
hold(ax1,'on')
ax2 = subplot(2,1,2);
hold(ax2,'on')

% use data from [gocue-0.5 to lastLick+0.5]


for i = 1:numel(dlc)
    cla(ax1)
    cla(ax2)

    tm = dlc(i).frameTimes;
    gc = gocueTimes(i);
    lastLick = lickL{i}(end) + 0.5;

    edges = [gc-0.2 lastLick+0.2];
    tix = findTimeIX(tm,edges,1);
    tm = tm(tix);
    
    % trialjaw = dlc(i).ts(:,coord,jawix);
    tjx = dlc(i).ts(:,1,jawix);
    tjy = dlc(i).ts(:,2,jawix);



    tjx = tjx(tix);
    tjy = tjy(tix);
    
    trialjaw = sqrt(tjx.^2 + tjy.^2);


    % fillnans
    trialjaw = fillmissing(trialjaw,'nearest');
    
    % Normalize
    trialjaw = (trialjaw - nanmean(trialjaw)) / nanstd(trialjaw);

    % hilbert transform works best on narrow band signals
    [b,a] = butter(3,50/400/2, 'low'); % [0.1/fs/2 40/fs/2] with islocalminima on phase works well at finding jaw closing times

    trialjaw = filtfilt(b,a,trialjaw);

    % Compute the analytical signal using the Hilbert transform
    analyticalSignal = hilbert(trialjaw);

    phase = angle(analyticalSignal);

    pvel = (gradient(phase));


    % find local minima in phase to get individual licks
    lickStartIX = findZeroCrossings(phase);
    lickStartIX = islocalmin(pvel,'FlatSelection','all','MaxNumExtrema',7);
    restingPosition = findRestingJawPosition(trialjaw);


    % plot for testing
    cm1 = colormap(flip(copper(ceil(numel(tm)/2))));
    cm2 = flip(cm1);
    cm = cat(1,cm1,cm2);
    colorData = (phase - min(phase)) / (max(phase) - min(phase));

    plot(ax1,tm,trialjaw,'k')
    s1=scatter(ax1,tm,trialjaw,20,colorData,'filled');
    s11=scatter(ax1,tm(lickStartIX),trialjaw(lickStartIX),40,'red','filled');
    colormap(cm)

    plot(ax2,tm,phase,'k');
    s2=scatter(ax2,tm,phase,20,colorData,'filled');
    s22=scatter(ax2,tm(lickStartIX),phase(lickStartIX),40,'red','filled');
    colormap(cm)

    xline(ax1,gocueTimes(i))
    xline(ax2,gocueTimes(i))

    % xlim(ax1,[gocueTimes(i)-0.2,ax1.XLim(2)])
    % xlim(ax2,[gocueTimes(i)-0.2,ax2.XLim(2)])

    xlabel(ax2,'Time from go cue (s)')
    ylabel(ax1,'Jaw Pos (px)')
    ylabel(ax2,'Phase (rad)')

    pause
end

%%

end