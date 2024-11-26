clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

% calculate CD choice for tagged units and single units, bootstrapped
% calculate CD choice selectivity matrices and transition latency b/w delay
% and response epochs
% see economo 2018, ED fig 7

%% LOAD DATA

datapth = 'C:\Users\munib\Documents\Economo-Lab\data\curated\';

alignEvent = 'goCue'; % 'goCue' 'lastLick' 'firstLick'

ptnfn = findMostRecentFile(fullfile(datapth,'PTNs'), alignEvent);
unfn = findMostRecentFile(fullfile(datapth,'UnidentifiedNeurons'), alignEvent);

ptn = load(ptnfn);
un = load(unfn);

% ptn and un should have save number of entries, each corresponds to a session

%%

clearvars -except alignEvent datapth ptnfn unfn utilspth ptn un

%% BOOTSTRAP PARAMS

boot.nBoots = 100;
boot.condLabels = {'rhit','lhit'};
boot.nSessions = 5;
boot.nUnits = 50;
boot.nTrials = 50;
% boot.cd.ix = findTimeIX(un.sessobj(1).time,[un.sesspar(1).eventTimes.sample un.sesspar(1).eventTimes.sample+0.5],1); % stimulus
% boot.cd.ix = findTimeIX(un.sessobj(1).time,[un.sesspar(1).eventTimes.goCue-0.6 un.sesspar(1).eventTimes.goCue],1); % choice
boot.cd.ix = findTimeIX(un.sessobj(1).time,[un.sesspar(1).eventTimes.goCue un.sesspar(1).eventTimes.goCue+0.5],1); % action
boot.pref_dir = 1;


nSessions = numel(ptn.meta);
nUnitsPerSession = floor(boot.nUnits/boot.nSessions);
nTime = numel(un.sessobj(1).time);

% baselineix = findTimeIX(un.sessobj(1).time,[-2.145,-1.85]);
baselineix = findTimeIX(un.sessobj(1).time,[-2 2]);

nTagged = 0;
for isess = 1:nSessions
    nTagged = nTagged + sum(ptn.tag(isess).nTag);
end

%% CODING DIRECTION (CHOICE)

%% TAGGED UNITS

clear tag %proj cd
tag.trialdat = cell(nTagged,1);
tag.session = zeros(nTagged,1);
% tag.psth = zeros(515,53,2); % (time,units,cond)
ct = 1;
for isess = 1:nSessions
    temptag = ptn.tag(isess);
    condix = find(ismember(ptn.sesspar(isess).condLabel,boot.condLabels));
    trix = balanceAndSplitTrials(ptn.sesspar(isess).trialid, condix, 1, 0); % {rhit,lhit}
    % trix = ptn.sesspar(isess).trialid(condix);
    for iprobe = 1:numel(temptag.nTag)
        for iunit = 1:temptag.nTag(iprobe)
            tag.trialdat{ct} = squeeze(temptag.trialdat{iprobe}(:,iunit,:));
            tag.session(ct) = isess;
            % for icond = 1:numel(trix)
            %     temp = nanmean(tag.trialdat{ct}(:,trix{icond}),2);
            %     tag.psth(:,ct,icond) = temp;
            % end
            ct = ct + 1;
        end
    end
end

% % mean subtract PSTH using baseline mean
% mu_ptn = nanmean(tag.psth(baselineix,:,:),[1 3]);
% mu_ptn = fillmissing(mu_ptn,"constant",nanmedian(mu_ptn));
% mu_ptn(mu_ptn==0) = 1;
% sigma_ptn = nanstd(tag.psth(baselineix,:,:),[],[1 3]);
% sigma_ptn = fillmissing(sigma_ptn,"constant",nanmedian(sigma_ptn));
% sigma_ptn(sigma_ptn==0) = 1;
%
% tag.psth = (tag.psth - mu_ptn) ./ sigma_ptn;
%
% if boot.pref_dir
%     % reorder psth so that preferred condition is first
%     a = tag.psth(boot.cd.ix,:,1) > tag.psth(boot.cd.ix,:,2);
%     a = a(:); % num time points for calculating cd in cond1 that are greater than that of cond2
%     a = sum(a) / numel(a);
%     if a < 0.5
%         temp = tag.psth;
%         tag.psth(:,:,1) = temp(:,:,2);
%         tag.psth(:,:,2) = temp(:,:,1);
%     end
% end
%
% % CALCULATE CD CHOICE
% cd.ptn = calcCD(tag.psth,boot.cd.ix);
% proj.ptn = tensorprod(tag.psth,cd.ptn,2,1);
%
% % figure; plot(mySmooth(proj.ptn,21,'reflect'))

%% SINGLE UNITS (BOOTSTRAP)

clear proj

proj.un = zeros(nTime,2,boot.nBoots); % (time,cond [rhit,lhit],bootstraps)
proj.ptn = zeros(nTime,2,boot.nBoots);
for iboot = 1:boot.nBoots
    clear psth

    if mod(iboot,20)==0; disp(['Iteration ' num2str(iboot) '/' num2str(boot.nBoots)]); end

    % GATHER SINGLE AND TAGGED UNITS

    % preallocate psth
    psth.un = nan(nTime,boot.nUnits,2);
    psth.ptn = nan(nTime,boot.nUnits,2);

    % sessions to sample
    sess2sample = randsample(nSessions,boot.nSessions,true);

    % for each session
    for isess = 1:boot.nSessions
        sessobj = un.sessobj(sess2sample(isess));
        sesspar = un.sesspar(sess2sample(isess));

        % get balanced number of trials
        condix = find(ismember(sesspar.condLabel,boot.condLabels));
        trials = balanceAndSplitTrials(sesspar.trialid, condix, 1, 0); % {rhit,lhit}
        trials = cellfun(@(x)  randsample(x,boot.nTrials,true), trials, 'uni', 0);
        % trials = sesspar.trialid(condix);

        % sample single units
        trialdat_un = cat(2,sessobj.trialdat{:});
        unitix = randsample(size(trialdat_un,2),nUnitsPerSession,true);
        trialdat_un = trialdat_un(:,unitix,:);

        % calculate PSTH
        for icond = 1:numel(trials)
            d = (1:nUnitsPerSession) + (nUnitsPerSession*(isess-1));
            psth.un(:,d,icond) = nanmean(trialdat_un(:,:,trials{icond}),3);
        end

        % mean subtract PSTH using baseline mean
        mu_un = nanmean(psth.un(baselineix,:,:),[1 3]);
        mu_un = fillmissing(mu_un,"constant",nanmedian(mu_un));
        mu_un(mu_un==0) = 1;
        sigma_un = nanstd(psth.un(baselineix,:,:),[],[1 3]);
        sigma_un = fillmissing(sigma_un,"constant",nanmedian(sigma_un));
        sigma_un(sigma_un==0) = 1;
    end
    psth.un = (psth.un - mu_un) ./ sigma_un;
    if boot.pref_dir
        % reorder psth so that preferred condition is first
        a = psth.un(boot.cd.ix,:,1) > psth.un(boot.cd.ix,:,2);
        a = a(:); % num time points for calculating cd in cond1 that are greater than that of cond2
        a = sum(a) / numel(a);
        if a < 0.5
            temp = psth.un;
            psth.un(:,:,1) = temp(:,:,2);
            psth.un(:,:,2) = temp(:,:,1);
        end
    end

    % GATHER TAGGED UNITS

    % tagged units to sample
    tag2sample = randsample(nTagged,boot.nUnits,true);
    boottag.trialdat = tag.trialdat(tag2sample);
    boottag.session = tag.session(tag2sample);

    % for each session
    for iunit = 1:boot.nUnits
        % sample tagged units
        trialdat_ptn = boottag.trialdat{iunit};
        thissess = boottag.session(iunit);
        sesspar = ptn.sesspar(thissess);

        % get balanced number of trials
        condix = find(ismember(sesspar.condLabel,boot.condLabels));
        trials = balanceAndSplitTrials(sesspar.trialid, condix, 1, 0); % {rhit,lhit}
        trials = cellfun(@(x)  randsample(x,boot.nTrials,true), trials, 'uni', 0);

        % calculate PSTH
        for icond = 1:numel(trials)
            psth.ptn(:,iunit,icond) = nanmean(trialdat_ptn(:,trials{icond}),2);
        end
    end

    % mean subtract PSTH using baseline mean
    mu_ptn = nanmean(psth.ptn(baselineix,:,:),[1 3]);
    mu_ptn = fillmissing(mu_ptn,"constant",nanmedian(mu_ptn));
    sigma_ptn = nanstd(psth.ptn(baselineix,:,:),[],[1 3]);
    sigma_ptn = fillmissing(sigma_ptn,"constant",nanmedian(sigma_ptn));

    psth.ptn = (psth.ptn - mu_ptn) ./ sigma_ptn;

    if boot.pref_dir
        % reorder psth so that preferred condition is first
        a = psth.ptn(boot.cd.ix,:,1) > psth.ptn(boot.cd.ix,:,2);
        a = a(:); % num time points for calculating cd in cond1 that are greater than that of cond2
        a = sum(a) / numel(a);
        if a < 0.5
            temp = psth.ptn;
            psth.ptn(:,:,1) = temp(:,:,2);
            psth.ptn(:,:,2) = temp(:,:,1);
        end
    end

    % CALCULATE CD CHOICE
    cd.un = calcCD(psth.un,boot.cd.ix);
    cd.ptn = calcCD(psth.ptn,boot.cd.ix);

    proj.un(:,:,iboot) = tensorprod(psth.un,cd.un,2,1);
    proj.ptn(:,:,iboot) = tensorprod(psth.ptn,cd.ptn,2,1);

    % figure; plot(proj.un(:,:,iboot)); hold on; plot(proj.ptn(:,:,iboot))

end

%% PLOT CD SELECTIVITY FOR TAGGED AND SINGLE UNITS
clear sel

sel.un.data = abs(squeeze(proj.un(:,1,:) - proj.un(:,2,:)));
sel.ptn.data = abs(squeeze(proj.ptn(:,1,:) - proj.ptn(:,2,:)));

[sel.un.mu,sel.un.ci95] = mean_CI(sel.un.data,0.95,true);
[sel.ptn.mu,sel.ptn.ci95] = mean_CI(sel.ptn.data,0.95,true);

sel.un.mu = sel.un.mu - mean(sel.un.mu(baselineix,:));
sel.ptn.mu = sel.ptn.mu - mean(sel.ptn.mu(baselineix,:));

cols = getColors;
f = figure;
f.Position = [300   272   570   343];
ax = prettifyAxis(gca);
hold on;
shadedErrorBar(un.sessobj(1).time,sel.un.mu,sel.un.ci95./1.5,{'Color',cols.un,'LineWidth',2},0.3,ax);
shadedErrorBar(un.sessobj(1).time,sel.ptn.mu,sel.ptn.ci95./2,{'Color',cols.ptn,'LineWidth',2},0.3,ax);
% shadedErrorBar(un.sessobj(1).time,sel.un.mu,sel.un.ci95,{'Color',cols.un,'LineWidth',2},0.3,ax);
% shadedErrorBar(un.sessobj(1).time,sel.ptn.mu,sel.ptn.ci95,{'Color',cols.ptn,'LineWidth',2},0.3,ax);
xlabel('Time from go cue (s)')
ylabel('Sel. in CD (a.u.)')
xlim([-2.15 2.5])
plotEventTimes(ax,un.sesspar(1).eventTimes)
plot(ax.XLim,[0 0],'k-')


leg = legend('Single units', 'PTNs'); leg.Box = 'off';




































