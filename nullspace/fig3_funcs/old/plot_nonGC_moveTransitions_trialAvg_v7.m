function plot_nonGC_moveTransitions_trialAvg_v7(dat,obj,me,rez,params,meta,ndims)

rng(pi)


% m2q
for sessix = 1:numel(me)
    f = figure;
%     f.Position = [403   412   809   543];
    ax = gca;
    hold on

    dt = params(sessix).dt;

    % transition data
    qdat = dat.qdat{sessix};
    nTrials = numel(qdat);


    % sort null/potent spaces by ve
    [~,nullix] = sort(rez(sessix).ve.null,'descend');
    nullix = nullix(1:ndims);

    [~,potentix] = sort(rez(sessix).ve.potent,'descend');
    potentix = potentix(1:ndims);

    nbouts = 5; % approximately max number of bouts in a trial
    allnull = nan(500,nTrials*nbouts,ndims);
    allpotent = nan(500,nTrials*nbouts,ndims);
    allme = nan(500,nTrials*nbouts,ndims);

    for dimix = 1:ndims

        for trix = 1:nTrials
            qdat_trix = qdat{trix};
            if isempty(qdat_trix)
                continue
            end

            sm = 41;
            null = mySmooth(rez(sessix).N_null(:,trix,nullix(dimix)),sm);
            potent = mySmooth(rez(sessix).N_potent(:,trix,potentix(dimix)),sm);
            move = mySmooth(me(sessix).move(:,trix),sm);
            medat = me(sessix).data(:,trix);

            % for each transition
            for i = 1:size(qdat_trix,2)
                tix = qdat_trix(1,i):qdat_trix(end,i);
                ix = find(tix == 0);
                if ~isempty(ix)
                    tix(ix) = [];
                end
                ts = tix .* dt;
                ts = ts - mean(ts);

                nullplot = null(tix);
                potentplot = potent(tix);
                meplot = medat(tix);

                % get center of tix
                center = ceil(numel(tix)/2);
                try
                    allnull(250-center:250+center-1,(trix+i-1),dimix) = nullplot;
                    allpotent(250-center:250+center-1,(trix+i)-1,dimix) = potentplot;
                    allme(250-center:250+center-1,(trix+i)-1,dimix) = meplot;
                catch
                    allnull(250-center:250+center-2,(trix+i-1),dimix) = nullplot;
                    allpotent(250-center:250+center-2,(trix+i)-1,dimix) = potentplot;
                    allme(250-center:250+center-2,(trix+i)-1,dimix) = meplot;
                end


            end

        end


    end

     % sumsqmag
    meannull = sum(allnull.^2,3,'omitnan');
    meanpotent = sum(allpotent.^2,3,'omitnan');
%     meannull = squeeze(nanmean(allnull.^2,3));
%     meanpotent = squeeze(nanmean(allpotent.^2,3));


    nullmeanplot = nanmean(meannull,2);
    nullerrplot = nanstd(meannull,[],2) ./ sqrt(size(allnull,2));
    potentmeanplot = nanmean(meanpotent,2);
    potenterrplot = nanstd(meanpotent,[],2) ./ sqrt(size(allnull,2));

    yyaxis(ax,'left')
    alph = 0.1;
     shadedErrorBar((1:500)./100 - 2.5,nullmeanplot,nullerrplot,{'Color',[62, 168, 105]./255,'LineWidth',2,'LineStyle','-'},alph,ax)
    shadedErrorBar((1:500)./100 - 2.5,potentmeanplot,potenterrplot,{'Color',[255, 56, 140]./255,'LineWidth',2,'LineStyle','-'},alph,ax)
%     shadedErrorBar(obj(sessix).time,nullmeanplot,nullerrplot,{'Color',[62, 168, 105]./255,'LineWidth',2,'LineStyle','-'},0.5,ax)
%     shadedErrorBar(obj(sessix).time,potentmeanplot,potenterrplot,{'Color',[255, 56, 140]./255,'LineWidth',2,'LineStyle','-'},0.5,ax)
    ax.YColor = [100 100 100 ]./255;
    xlabel('Time to move (s)')
    ylabel('SSM - all dims (a.u)')

    ylims = ax.YLim;

    meanme = squeeze(normalize(nanmean(allme,2),'range',[ylims(1) ylims(2)]));
%     meanme = squeeze(nanmean(allme,2));
    meanme = meanme(:,1);
    yyaxis(ax,'right')
    plot((1:500)./100 - 2.5,meanme,'Color',[225, 144, 15]./255,'LineWidth',3);
    ax.YColor = [225, 144, 15]./255;
    ylabel('Motion Energy')

    xline(0,'k:','LineWidth',2)

    xlim([-0.25 0.25])
    title([meta(sessix).anm ' - ' meta(sessix).date])
    

    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'-','Color',[225, 144, 15]./255);
    h(2) = plot(NaN,NaN,'-','Color',[62, 168, 105]./255);
    h(3) = plot(NaN,NaN,'-','Color',[255, 56, 140]./255);
    legend(h, 'motion energy','null','potent');
    
%     break
end



end