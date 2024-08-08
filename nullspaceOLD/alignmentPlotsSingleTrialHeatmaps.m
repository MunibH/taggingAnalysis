%% NP alignment and PSTHs

% find L/R selective cells per session
edges = [-0.8 0];

evtimes = getEventTimes(obj(1).bp.ev,{'bitStart','sample','delay','goCue'},'goCue');

close all
for sessix = 1:numel(obj)
    cond2use = [5 6]; %[8 9];
    cluix = findSelectiveCells(obj(sessix),params(sessix),edges,cond2use);

    cond2use = {'hit|miss'};
    clear trialdat W proj r2
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix);

    trialdat.full = permute(obj(sessix).trialdat,[1 3 2]);
    trialdat.full = trialdat.full(:,trix,:);

    W.null = rez(sessix).Qnull;
    W.potent = rez(sessix).Qpotent;
    fns = {'null','potent'};
    psth = [];
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.recon.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)

        for k = 1:numel(clus) % for each cell
            thisclu = clus(k);
            if j==1
                psth = cat(2,psth,obj(sessix).psth(:,thisclu,[5,6]));
            end
            % calculcate variance explained by CD choice
            orig = trialdat.full(:,:,thisclu); % (time,trials)

            fr = mean(mean(orig)); % subspace contribution method
            % weight = norm(W.(fns{j})(k));
            % r2.(fns{j}){sessix}(k) = fr*weight;

            recon = trialdat.recon.(fns{j})(:,:,thisclu); % (time,trials) % ve by recon method
            mdl = fitlm(orig(:),recon(:));
            r2.(fns{j})(k) = mdl.Rsquared.Ordinary;
        end
    end
    alignment = (r2.null - r2.potent) ./ (r2.null + r2.potent);
    [alignment_,sortix] = sort(alignment);
    clus = clus(sortix);
    
    pmask = alignment_ <= -0.6;
    nmask = alignment_ >= 0.6;
    if sum(pmask)==0 || sum(nmask)==0
        continue
    end
    clus_ = clus(nmask | pmask);
    alignment_ = alignment_(nmask | pmask);
    data = permute(obj(sessix).trialdat(:,clus_,:),[1 3 2]);
    f = figure;
    t = tiledlayout('flow');

    cond2use = [5 6];
    trix = cell2mat(params(sessix).trialid(cond2use)');
    nTrials = numel(trix);
    data = data(:,trix,:);
    for i = 1:numel(alignment_)
        ax = nexttile;
        hold on;

        this = data(:,:,i);

        imagesc(obj(sessix).time,1:nTrials,this')
        title(['Alignment = ' num2str(round(alignment_(i),2))])
        xlabel('Time from go cue (s)')
        ylabel('Trials')
        plotEventTimes(ax,evtimes)
        xlim([-2.3,2])

    end
    

    % pdata = sum(data(:,:,pmask).^2,3);
    % ndata = sum(data(:,:,nmask).^2,3);

    % f = figure;
    % ax = subplot(1,2,1);
    % hold on;
    % imagesc(obj(1).time,1:size(pdata,2),pdata')
    % title('Potent-aligned units')
    % xlabel('Time from go cue (s)')
    % ylabel('Trials')
    % plotEventTimes(ax,evtimes)
    % xlim([-2.3,2])
    % 
    % ax = subplot(1,2,2);
    % hold on;
    % imagesc(obj(1).time,1:size(ndata,2),ndata')
    % title('Null-aligned units')
    % plotEventTimes(ax,evtimes)
    % xlim([-2.3,2])

    % sel = psth(:,:,1) - psth(:,:,2);
    % [~,sortix] = sort(alignment);
    % f = figure;
    % ax = gca;
    % hold on;
    % imagesc(obj(sessix).time,1:size(sel,2),sel')
    % % title(['Alignment = ' num2str(round(alignment(i),2))])
    % xlabel('Time from go cue (s)')
    % ylabel('Units')
    % plotEventTimes(ax,evtimes)
    % xlim([-2.3,2])

    % % PLOT PSTHS AND ALIGNMENT
    % f = figure;
    % t = tiledlayout('flow');
    % for i = 1:numel(alignment)
    %     ax = nexttile;
    %     hold on;
    %     ax = prettifyPlot(ax);
    % 
    %     thisclu = clus(i);
    % 
    %     plot(obj(sessix).time,obj(sessix).psth(:,thisclu,5),'b')
    %     plot(obj(sessix).time,obj(sessix).psth(:,thisclu,6),'r')
    %     title(['Alignment = ' num2str(round(alignment(i),2))])
    %     xlabel('Time from go cue (s)')
    %     ylabel('Firing rate')
    %     plotEventTimes(ax,evtimes)
    %     xlim([-2.3,2])
    % 
    % end

end


%% NP alignment and Depth

% make it work for 2 probe sessions
% get data for nan sessions

cols = linspecer(numel(meta));

f = figure;
ax = gca;
hold on;
for isess = 1:numel(meta)
    clear d
    % get depth
    depth = getDepth(meta(isess),obj(isess));
    % subset by clus of interest
    if ~iscell(depth)
        if isnan(depth)
            disp([meta(isess).anm ' ' meta(isess).date])
            continue
        end
    end
    if numel(meta(isess).probe) > 1
        for iprobe = 1:numel(meta(isess).probe)
            d{iprobe} = depth{iprobe}(params(isess).cluid{iprobe});
            plot(d{iprobe},'Color',cols(isess,:));
            plot(d{iprobe},'.','MarkerSize',10,'Color',cols(isess,:))
        end
        continue
    end

    d = depth(params(isess).cluid);
    plot(d,'Color',cols(isess,:));
    plot(d,'.','MarkerSize',10,'Color',cols(isess,:))

end
xlabel('Unit num')
ylabel('Depth (um)')

%% CD Choice NP alignment and PSTHs












