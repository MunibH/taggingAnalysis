clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% PARAMETERS

params = defaultParams();

% % specify changes here

params.region = 'any'; % 'alm','tjm1','mc', 'any'
params.probeType = 'any'; % 'h2','np2','np1', 'any'

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);

% meta = loadJPV8(meta,datapth);
meta = loadJPV11(meta,datapth);
% meta = loadJPV12(meta,datapth);
% meta = loadJPV13(meta,datapth);

% subset meta
meta = meta(1);
% meta = subsetMetaByParams(meta,params);

%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

me = loadMotionEnergy(obj, meta, params, datapth);

kin = getKinematics(obj, me, params);

% TAGGED UNIT META

tag.nTag = numel(obj.tag);
tag.cluid.clu = [obj.tag(:).cluid]; % where tagged units are in obj.clu
tag.cluid.obj = find(ismember(params.cluid,tag.cluid.clu))'; % where in params.cluid, trialdat, psth

%% SELECTIVITY

% doing this with PSTHs for now, need to do it with single trials (balanced)

%<XY>: X=trial type, Y=lick direction

% conds
cond.RR = 2;
cond.LL = 3;
cond.RL = 4;
cond.LR = 5;

% trials per condition
fnames = fieldnames(cond);
for i = 1:numel(fnames)
    f = fnames{i};
    trials.(f) = params.trialid{cond.(f)};
end

% stimulus selectivity = <RR>+<RL> - <LL>+<LR>
sel.stim = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.RL) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.LR) );

% choice selectivity = <RR>+<LR> - <LL>+<RL>
sel.choice = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LR) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.RL) );

% action selectivity = <RR>+<LR> - <LL>+<RL>
sel.action = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LR) ) - ...
    ( obj.psth(:,:,cond.LL)+obj.psth(:,:,cond.RL) );

% outcome selectivity = <LL>+<RR> - <LR>+<RL>
sel.outcome = ( obj.psth(:,:,cond.RR)+obj.psth(:,:,cond.LL) ) - ...
    ( obj.psth(:,:,cond.LR)+obj.psth(:,:,cond.RL) );

close all
f = figure;
ax = prettifyAxis(gca);
hold on;
fnames = fieldnames(sel);
cols = linspecer(numel(fnames));
for i = 1:numel(fnames)
    f = fnames{i};
    mu = mean(sel.(f),2);
    ci = getCI(sel.(f));
    shadedErrorBar(obj.time,mu,ci,{'Color',cols(i,:),'LineWidth',2,'HandleVisibility','off'},0.2,ax)
end
xlim([obj.time(1),obj.time(end)])
plotEventTimes(ax,params.eventTimes)


%% find selective neurons
clear times mu sel cd

%<XY>: X=trial type, Y=lick direction

nboots = 1;

% conds
cond.RR = 2;
cond.LL = 3;
cond.RL = 4;
cond.LR = 5;
% cond.RR_LL = 6; % just concatening sampled RR and LL so that they're using the same trials



% times (relative to go cue)
times.stimulus   = [-1.85 -1.2];
times.choice  = [-0.6 0];
times.action  = [0 0.6];
% times.outcome = [2 2.995];
times.presample = [-2.14 -1.85];
times.prego = [-0.1 0];
times.postgo = [0 0.1];
fnames = fieldnames(times);
for i = 1:numel(fnames)
    f = fnames{i};
    times.ix.(f) = findTimeIX(obj.time,times.(f), true);
end

% get avg baseline activity across all conditions
alltrials = cat(1,trials.RR,trials.LL,trials.RL,trials.LR);
baseline_data = squeeze(mean(obj.trialdat(times.ix.presample,:,alltrials), 1));
baseline_mean = mean(baseline_data,2)'; % baseline mean across trials (neurons,1)
% baseline_std = std(baseline_data,[],2)'; % baseline std across trials
baseline_std = std(squeeze(mean(obj.trialdat(:,:,alltrials), 1)),[],2)';
baseline_std(baseline_std==0) = 1;
trialdat_norm = (obj.trialdat - baseline_mean);% ./ baseline_std;

% calculate cds/selectivity/dprime
for iboot = 1:nboots
    if mod(iboot,50)==0
        disp(['Iteration: ' num2str(iboot) '/' num2str(nboots)])
    end
    % balance # of correct trials and balance # of error trials
    correct_trials = balanceAndSplitTrials(params.trialid,[cond.RR,cond.LL],0.7,0.3);
    trials.RR = correct_trials.train{1};
    trials.LL = correct_trials.train{2};

    error_trials = balanceAndSplitTrials(params.trialid,[cond.RL,cond.LR],0.9,0.1);
    trials.RL = error_trials.train{1};
    trials.LR = error_trials.train{2};
    
    trials.RR_LL = cat(1,trials.RR,trials.LL);

    epochs = fieldnames(times.ix);
    for iepoch = 1:numel(epochs) % stim,choice,action,outcome
        e = epochs{iepoch};
        condnames = fieldnames(trials);
        for icond = 1:numel(condnames)
            f = condnames{icond};
            thesetrials = trials.(f);
    
            N = numel(thesetrials);
    
            T = numel(times.ix.(e));
    
            fr = trialdat_norm(times.ix.(e),:,thesetrials);
            fr_sum_trials = squeeze(sum(fr,3));
            fr_sum_time_trials = sum(fr_sum_trials,1)';
    
            mu.(e).(f) = 1/(N*T) * fr_sum_time_trials; % (neurons,1)
        end
    end
    
    % % sel.stim(:,iboot,1) = mu.RR+mu.RL; % (neurons,boots,group)
    % % sel.stim(:,iboot,2) = mu.LL+mu.LR; % (neurons,boots,group)
    % 
    % sel.stim(:,iboot,1) = mu.RR; % (neurons,boots,group)
    % sel.stim(:,iboot,2) = mu.LL; % (neurons,boots,group)
    % 
    % cd(:,iboot) = (mu.RR+mu.RL) - (mu.LL+mu.LR);
    % % cd(:,iboot) = mu.RR - mu.LL;

    cd.stimulus(:,iboot) = (mu.stimulus.RR+mu.stimulus.RL) - (mu.stimulus.LL+mu.stimulus.LR);
    cd.choice(:,iboot) = (mu.choice.RR+mu.choice.LR) - (mu.choice.LL+mu.choice.RL);
    cd.action(:,iboot) = (mu.action.RR+mu.action.LR) - (mu.action.LL+mu.action.RL);
    % cd.outcome(:,iboot) = (mu.outcome.RR+mu.outcome.LL) - (mu.outcome.RL+mu.outcome.LR);
    cd.ramping(:,iboot) = mu.choice.RR_LL - mu.presample.RR_LL;
    % cd.ramping(:,iboot) = mu.postgo.RR_LL - mu.presample.RR_LL;
    cd.go(:,iboot) = mu.postgo.RR_LL - mu.prego.RR_LL;

end

% 
% %%
% for i = 1:numel(params.cluid)
%     % [p(i),h(i)] = ranksum(squeeze(sel.stim(i,:,1)) , squeeze(sel.stim(i,:,2)) , 'alpha',0.01 );
%     [h(i),p(i)] = ttest(squeeze(sel.stim(i,:,1)) , squeeze(sel.stim(i,:,2)) , 'alpha',0.001);
% end
% h(isnan(h)) = 0;
% find(~h)
% 
% 
% 
% %%
% a = squeeze(sel.stim(9,:,:));
% figure; plot(a)


%% project single trials onto CD
clear proj recon
cdnames = fieldnames(cd);

% convert cd to matrix
cdmat = [];
for icd = 1:numel(cdnames)
    cdmat = cat(2,cdmat,cd.(cdnames{icd}));
end
% orthogonalize
cdmat = gschmidt(cdmat);

% proj
proj = tensorprod(trialdat_norm,cdmat,2,1);

% recon and ve (for each neuron and cd)
dat = permute(trialdat_norm,[1 3 2]); % (time,trials,neurons)
for icd = 1:numel(cdnames)
    thiscd = cdnames{icd};
    recon.(thiscd) = tensorprod(proj(:,:,icd),cdmat(:,icd),3,2);
    for iunit = 1:size(dat,3)
        temp1 = mean(dat(:,:,iunit),2);
        temp2 = mean(recon.(thiscd)(:,:,iunit),2);
        cc = corrcoef(temp1(:),temp2(:));
        r2(iunit,icd) = cc(1,2).^2;
    end
end


%% plot cd proj

close all
f = figure;
t = tiledlayout('flow');
for icd = 1:numel(cdnames)
    ax = prettifyAxis(nexttile);
    hold on;
    thisproj = proj(:,:,icd); % (time,trials)
    r = thisproj(:,params.trialid{2});
    l = thisproj(:,params.trialid{3});

    mur = mean(r,2);
    mul = mean(l,2);
    er = getCI(r);
    el = getCI(l);
    shadedErrorBar(obj.time,mur,er,{'Color','b','LineWidth',1},0.2,ax);
    shadedErrorBar(obj.time,mul,el,{'Color','r','LineWidth',1},0.2,ax);
    title(ax,cdnames{icd},'fontsize',10.5)
    plotEventTimes(ax,params.eventTimes)
end
xlabel(t,'Time from go cue (s)')
ylabel(t,'Proj (a.u.)')


%% plot unit and reconstruction

for icd = 1:numel(cdnames)
    thiscd = cdnames{icd};
    recon_plus_mean.(thiscd) = permute(recon.(thiscd),[1 3 2]) + baseline_mean;
end

f = figure;
ax = prettifyAxis(gca);
hold on;
for i = 1:numel(params.cluid)
    cla(ax)

    thisdat = squeeze(obj.trialdat(:,i,:));
    
    r = thisdat(:,params.trialid{2});
    l = thisdat(:,params.trialid{3});

    mur = mean(r,2);
    mul = mean(l,2);
    er = getCI(r);
    el = getCI(l);
    shadedErrorBar(obj.time,mur,er,{'Color','b','LineWidth',1},1,ax);
    shadedErrorBar(obj.time,mul,el,{'Color','r','LineWidth',1},1,ax);

    thisrecon = squeeze(recon_plus_mean.action(:,i,:));

    r = thisrecon(:,params.trialid{2});
    l = thisrecon(:,params.trialid{3});

    mur = mean(r,2);
    mul = mean(l,2);
    er = getCI(r);
    el = getCI(l);
    shadedErrorBar(obj.time,mur,er,{'Color',[0,0.7,1],'LineWidth',1},1,ax);
    shadedErrorBar(obj.time,mul,el,{'Color',[1,0.2,0.7],'LineWidth',1},1,ax);

    plotEventTimes(ax,params.eventTimes,'k',false)

    pause
    
end



%% variance explained

X = reshape(trialdat_norm,size(trialdat_norm,1)*size(trialdat_norm,2),size(trialdat_norm,3));
muX = mean(X,1);
total_var_X = sum(var(X,[],1));

for icd = 1:numel(cdnames)
    thisrecon = reshape(recon.(cdnames{icd}),size(trialdat_norm,1)*size(trialdat_norm,2),size(trialdat_norm,3));
    residual_var = sum(var(X-thisrecon,[],1));
    ve.(cdnames{icd}) = (total_var_X - residual_var) / total_var_X;

end
%%
% null
Y = reshape(recon.null,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.null = (total_var_X - residual_var) / total_var_X;

% potent
Y = reshape(recon.potent,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.potent = (total_var_X - residual_var) / total_var_X;

% cd
Y = reshape(recon.cd,size(trialdat_zscored,1)*size(trialdat_zscored,2),size(trialdat_zscored,3));
residual_var = sum(var(X-Y,[],1));
ve.cd = (total_var_X - residual_var) / total_var_X;











