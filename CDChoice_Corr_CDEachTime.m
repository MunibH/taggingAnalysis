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
boot.cd.ix = findTimeIX(un.sessobj(1).time,[un.sesspar(1).eventTimes.goCue-0.6 un.sesspar(1).eventTimes.goCue],1); % choice
% boot.cd.ix = findTimeIX(un.sessobj(1).time,[un.sesspar(1).eventTimes.goCue un.sesspar(1).eventTimes.goCue+0.5],1); % action
boot.pref_dir = 1;


nSessions = numel(ptn.meta);
nUnitsPerSession = floor(boot.nUnits/boot.nSessions);
nTime = numel(un.sessobj(1).time);

baselineix = findTimeIX(un.sessobj(1).time,[-2.145,-1.85]);
% baselineix = findTimeIX(un.sessobj(1).time,[-2 2]);

nTagged = 0;
for isess = 1:nSessions
    nTagged = nTagged + sum(ptn.tag(isess).nTag);
end

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

%% SINGLE UNITS (BOOTSTRAP)

clear cd

normix = findTimeIX(ptn.tag(1).time,[-0.4,0],1);

for iboot = 1:boot.nBoots
    clear psth
    % disp(['Iteration ' num2str(iboot) '/' num2str(boot.nBoots)]);
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
    end
    % mean subtract PSTH using baseline mean
    mu_un = nanmean(psth.un(baselineix,:,:),[1 3]);
    mu_un = fillmissing(mu_un,"constant",nanmedian(mu_un));
    mu_un(mu_un==0) = 1;
    sigma_un = nanstd(psth.un(baselineix,:,:),[],[1 3]);
    sigma_un = fillmissing(sigma_un,"constant",nanmedian(sigma_un));
    sigma_un(sigma_un==0) = 1;

    psth.un = (psth.un - mu_un) ./ sigma_un;

    clear pref_dir temp
    if boot.pref_dir
        for iunit = 1:size(psth.un,2)
            temp = squeeze(psth.un(boot.cd.ix,iunit,:));
            a = temp(:,1) > temp(:,2);
            a = a(:);
            a = sum(a) / numel(a);
            if a < 0.5
                pref_dir(iunit,:) = [2 1];
            else
                pref_dir(iunit,:) = [1 2];
            end
        end
    end
    temp = psth.un;
    for iunit = 1:size(psth.un,2)
        temp(:,iunit,1) = psth.un(:,iunit,pref_dir(iunit,1));
        temp(:,iunit,2) = psth.un(:,iunit,pref_dir(iunit,2));
    end
    psth.un = temp;


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

    clear pref_dir temp
    if boot.pref_dir
        for iunit = 1:size(psth.ptn,2)
            temp = squeeze(psth.ptn(boot.cd.ix,iunit,:));
            a = temp(:,1) > temp(:,2);
            a = a(:);
            a = sum(a) / numel(a);
            if a < 0.5
                pref_dir(iunit,:) = [2 1];
            else
                pref_dir(iunit,:) = [1 2];
            end
        end
    end
    temp = psth.ptn;
    for iunit = 1:size(psth.ptn,2)
        temp(:,iunit,1) = psth.ptn(:,iunit,pref_dir(iunit,1));
        temp(:,iunit,2) = psth.ptn(:,iunit,pref_dir(iunit,2));
    end
    psth.ptn = temp;

    % CALCULATE CD CHOICE
    cd.un = calcCD(psth.un,boot.cd.ix);
    cd.ptn = calcCD(psth.ptn,boot.cd.ix);


    % CALCULATE CD EACH TIME POINT

    for itime = 1:size(psth.un,1)    
        cd.un_time(itime,:) = calcCD(psth.un(itime,:,:),1);
        cd.ptn_time(itime,:) = calcCD(psth.ptn(itime,:,:),1);
    end


    for itime = 1:size(psth.un,1) 
        cd.un_corr(itime,iboot) = corr(cd.un_time(itime,:)',cd.un);
        cd.ptn_corr(itime,iboot) = corr(cd.ptn_time(itime,:)',cd.ptn);
    end
    
    % cd.un_corr(:,iboot) = cd.un_corr(:,iboot) ./ nanmean(cd.un_corr(normix,iboot));
    % cd.ptn_corr(:,iboot) = cd.ptn_corr(:,iboot) ./ nanmean(cd.ptn_corr(normix,iboot));

end


%% PLOT

cols = getColors;

f = figure;
f.Position = [680   564   512   314];
ax = prettifyAxis(gca);
hold on;

[m,h] = mean_CI(cd.un_corr,0.5,true);
shadedErrorBar(ptn.tag(1).time*1000,m,h,{'Color',cols.un,'LineWidth',2},0.3,ax)
[m,h] = mean_CI(cd.ptn_corr,0.5,true);
shadedErrorBar(ptn.tag(1).time*1000,m,h,{'Color',cols.ptn,'LineWidth',2},0.3,ax)


xlim([-80 100])
xlabel('Time from go cue (ms)')
ylabel({'Normalized correlation','with CD choice'})
ylim([0 ax.YLim(2)])
plot([0 0],ax.YLim,'k--')



%% OLD

%% TAGGED UNITS

clear tag %proj cd
tag.trialdat = cell(nTagged,1);
tag.session = zeros(nTagged,1);
tag.psth = zeros(515,53,2); % (time,units,cond)
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
            for icond = 1:numel(trix)
                temp = nanmean(tag.trialdat{ct}(:,trix{icond}),2);
                tag.psth(:,ct,icond) = temp;
            end
            ct = ct + 1;
        end
    end
end

% mean subtract PSTH using baseline mean
mu_ptn = nanmean(tag.psth(baselineix,:,:),[1 3]);
mu_ptn = fillmissing(mu_ptn,"constant",nanmedian(mu_ptn));
mu_ptn(mu_ptn==0) = 1;
sigma_ptn = nanstd(tag.psth(baselineix,:,:),[],[1 3]);
sigma_ptn = fillmissing(sigma_ptn,"constant",nanmedian(sigma_ptn));
sigma_ptn(sigma_ptn==0) = 1;

tag.psth = (tag.psth - mu_ptn) ./ sigma_ptn;

if boot.pref_dir
    for iunit = 1:size(tag.psth,2)
        temp = squeeze(tag.psth(boot.cd.ix,iunit,:));
        a = temp(:,1) > temp(:,2);
        a = a(:);
        a = sum(a) / numel(a);
        if a < 0.5
            pref_dir(iunit,:) = [2 1];
        else
            pref_dir(iunit,:) = [1 2];
        end
    end
end
temp = tag.psth;
for iunit = 1:size(tag.psth,2)
    temp(:,iunit,1) = tag.psth(:,iunit,pref_dir(iunit,1));
    temp(:,iunit,2) = tag.psth(:,iunit,pref_dir(iunit,2));
end
tag.psth = temp;

%%

[nTime,nUnits,nCond] = size(tag.psth);

for itime = 1:nTime

    temp = tag.psth(itime,:,:);

    cd(itime,:) = calcCD(temp,1);


end

cdchoice = calcCD(tag.psth,boot.cd.ix);

%%


correlationMatrix = zeros(nTime, nTime);

% Loop over each pair of time points to calculate correlations
for t1 = 1:nTime
    for t2 = 1:nTime
        % Compute the correlation between cd(t1,:) and cd(t2,:)
        correlationMatrix(t1, t2) = corr(cd(t1, :)', cd(t2, :)');
    end
end


%%
cols = getColors;


for t1 = 1:nTime
    cc(t1) = corr(cd(t1,:)',cdchoice);
end

normix = findTimeIX(ptn.tag(1).time,[-0.4,0],1);
cc = cc ./ nanmean(cc(normix));

%%

f = figure;
f.Position = [680   564   512   314];
ax = prettifyAxis(gca);
hold on;
plot(ptn.tag(1).time*1000,cc,'Color',cols.ptn,'LineWidth',3)
xlim([-80 100])
xlabel('Time from go cue (ms)')
ylabel({'Normalized correlation','with CD choice'})
ylim([0 ax.YLim(2)])
plot([0 0],ax.YLim,'k--')








