%% make susu fig 3 plots

close all

cols = getColors();

cond2plot = [2,3];
trix = cell2mat(params.trialid(cond2plot)');
for i = 1:numel(cond2plot)
    ctrix{i} = params.trialid{cond2plot(i)};
end

f = figure;
f.Renderer = 'painters';
ax = prettifyPlot(gca);
hold on;
histogram(YY.R2,20,'Normalization','count')
xlabel('Test-R2')
ylabel('Unit count')


baselineedges = [params.eventTimes.bitStart+0.1 params.eventTimes.sample];
baselineix = findTimeIX(obj.time,baselineedges);
bix = baselineix(1):baselineix(2);


f = figure;
f.Renderer = 'painters';
ax1 = subplot(2,2,1);
hold on;
ax2 = subplot(2,2,2);
hold on;
ax3 = subplot(2,2,3);
hold on;
ax4 = subplot(2,2,4);
hold on;
for i = size(y,3)-1:size(y,3)%1:size(y,3)
    cla(ax1)
    cla(ax2)
    cla(ax3)
    cla(ax4)

    % if YY.R2(i)<0.3%sel.delay.cluid(i)
    %     continue
    % end

    tempy = y(:,:,i);
    tempyhat = yhat(:,:,i);

    % put tempy in same scale as tempyhat
    yhatstd = std(tempyhat,[],1);
    tempyhat = tempyhat ./ yhatstd;
    tempy = tempy;% .* yhatstd;

    % make sure the baselines start at 0
    basey = mean(tempy(bix,:),1);
    tempy = tempy - basey;
    baseyhat = mean(tempyhat(bix,:),1);
    tempyhat = tempyhat - baseyhat;


    hold(ax1,'on')
    imagesc(ax1,obj.time,1:size(y,2),tempy');
    ax1.YDir = "reverse";
    title(ax1,'data')
    plotEventTimes(ax1,params.eventTimes)
    ylim(ax1,[1 size(y,2)])
    colorbar(ax1);

    hold(ax2,'on')
    imagesc(ax2,obj.time,1:size(y,2),tempyhat');
    ax2.YDir = "reverse";
    title(ax2,'pred')
    plotEventTimes(ax2,params.eventTimes)
    ylim(ax2,[1 size(y,2)])
    colorbar(ax2);

    
    tempr = tempy(:,ctrix{1});
    templ = tempy(:,ctrix{2});
    temprhat = tempyhat(:,ctrix{1});
    templhat = tempyhat(:,ctrix{2});

    mu.y{1} = mean(tempr,2);
    mu.y{2} = mean(templ,2);
    mu.yhat{1} = mean(temprhat,2);
    mu.yhat{2} = mean(templhat,2);

    ci.y{1} = getCI(tempr);
    ci.y{2} = getCI(templ);
    ci.yhat{1} = getCI(tempr);
    ci.yhat{2} = getCI(templ);

    hold(ax3,'on')
    shadedErrorBar(obj.time,mu.y{1},ci.y{1},{'Color',cols.rhit,'LineWidth',2},alph,ax3)
    shadedErrorBar(obj.time,mu.y{2},ci.y{2},{'Color',cols.lhit,'LineWidth',2},alph,ax3)
    hold(ax4,'on')
    shadedErrorBar(obj.time,mu.yhat{1},ci.yhat{1},{'Color',cols.rhit_aw,'LineWidth',2},alph,ax4)
    shadedErrorBar(obj.time,mu.yhat{2},ci.yhat{2},{'Color',cols.lhit_aw,'LineWidth',2},alph,ax4)
    title(ax4,['R2=' num2str(round(YY.R2(i),2))])
    

    % ylims(1) = min(ax3.YLim(1),ax4.YLim(1));
    % ylims(2) = max(ax3.YLim(2),ax4.YLim(2));
    % % ylims(1) = ylims(1)*1.1;
    % % ylims(2) = ylims(2)*1.1;
    % 
    % ax3.YLim = ylims;
    % ax4.YLim = ylims;

    plotEventTimes(ax3,params.eventTimes,'k',false)
    plotEventTimes(ax4,params.eventTimes,'k',false)
    xlabel(ax3,'Time from go cue (s)')
    ylabel(ax3,'zscored FR')
    xlabel(ax4,'Time from go cue (s)')
    ylabel(ax4,'zscored FR')
    title(ax3,[meta.anm ' ' meta.date ' Unit ' num2str(i)])

    pause

    clear mu ci
end

% clear mu ci