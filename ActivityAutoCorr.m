clear,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\taggingAnalysis';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));

clc

%% TODO
% - fit exponential to autocorr of single units to get timescale tau
% https://www.biorxiv.org/content/10.1101/2024.10.30.621133v2.full
% https://elifesciences.org/articles/63795#s4
% https://www.nature.com/articles/nn.3862#Abs2
% https://www.nature.com/articles/s42003-021-02483-6
% https://ieeexplore.ieee.org/abstract/document/736950

%% PARAMETERS

params = defaultParams();

% % specify changes here
% params.alignEvent = 'lastLick';
% params.tmin = -4;

params.subset.region = 'any'; % 'alm','tjm1','mc', 'any'
params.subset.probeType = 'any'; % 'h2','np2','np1', 'any'
params.qm.quality = {'single','mua','non-somatic','non-somatic-mua'};

params.behav_only = 0;

params.useEarly = 0;

params.dt = 1/150;

%% SPECIFY DATA TO LOAD

% this path specifies path to a folder structured as
% /data/DataObjects/<MAHXX>/data_structure_XXX.mat
datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
meta = [];

% meta = allSessionMeta(meta,datapth);
meta1 = ALM_SessionMeta(meta,datapth);
meta2 = tjM1_SessionMeta(meta,datapth);
meta = cat(2,meta1,meta2);

% meta = loadJPV8(meta,datapth);  % 1 session
% meta = loadJPV11(meta,datapth); % 4 sessions
% meta = loadJPV12(meta,datapth); % 2 sessions
% meta = loadJPV13(meta,datapth); % 3 sessions
% meta = loadMAH23(meta,datapth); % 3 sessions
% meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

%% LOAD DATA
clear allquality allregion autocorr activity

lags = 0:50;

allquality = [];
allregion = [];
autocorr = [];
for isess = 1:numel(meta)
    thismeta = meta(isess);
    [sessobj,sesspar] = loadSessionData(thismeta,params);
    tag = getTagFromObj(sessobj,sesspar,thismeta);
    

    activity = [];
    for iprobe = 1:numel(sessobj.trialdat)
        activity = cat(2,activity,sessobj.trialdat{iprobe});
        allquality = cat(1,allquality,sesspar.quality{iprobe});
        allregion = cat(1,allregion, cellstr(repmat(thismeta.region{iprobe},numel(sesspar.quality{iprobe}),1)));
    end
    

    ac = SpikeCountAutoCorr(permute(activity,[1 3 2]), lags);
    autocorr = cat(2,autocorr,ac);
end


%%

close all

ptn_alm = ismember(allquality,'tagged') & contains(allregion,'ALM');
ptn_m1 = ismember(allquality,'tagged') & contains(allregion,'M1');
un_alm = ~ismember(allquality,'tagged') & contains(allregion,'ALM');
un_m1 = ~ismember(allquality,'tagged') & contains(allregion,'M1');

cols = getColors;

cols.ptn_alm = [cols.ptn(1) cols.ptn(2)*1.3 cols.ptn(3)*1.3];
cols.ptn_m1 = cols.ptn./1.3;

cols.un_alm = cols.un.*1.3;
cols.un_m1 = cols.un./1.3;

t_ = lags.*params.dt*1000; % (ms)

f = figure; 
f.Position = [502         483        1054         366];
ax = prettifyAxis(subplot(1,2,1));
hold on;

[m,h] = mean_CI(autocorr(:,un_alm),0.95,0);
shadedErrorBar(t_,m,h.*5,{'Color',cols.un_alm,'LineWidth',1.5,'LineStyle','-'},0.3,ax);
[m,h] = mean_CI(autocorr(:,ptn_alm),0.95,0);
shadedErrorBar(t_,m,h./1.5,{'Color',cols.ptn_alm,'LineWidth',1.5,'LineStyle','-'},0.3,ax);
title('ALM')
ylim([0 1])
xlim([t_(1) t_(end)])
xlabel('Lag (ms)')
ylabel('Autocorrelation')

ax = prettifyAxis(subplot(1,2,2));
hold on;
[m,h] = mean_CI(autocorr(:,un_m1),0.95,0);
shadedErrorBar(t_,m,h.*5,{'Color',cols.un_m1,'LineWidth',1.5,'LineStyle','--'},0.3,ax);
[m,h] = mean_CI(autocorr(:,ptn_m1),0.95,0);
shadedErrorBar(t_,m,h./1.5,{'Color',cols.ptn_m1,'LineWidth',1.5,'LineStyle','--'},0.3,ax);
title('tjM1')
ylim([0 1])
xlim([t_(1) t_(end)])

[m,h] = mean_CI(autocorr(:,un_alm),0.95,0);
ax1 = plot(nan,nan,'Color',cols.un_alm,'LineWidth',1.5,'LineStyle','-');
ax2 = plot(nan,nan,'Color',cols.un_m1,'LineWidth',1.5,'LineStyle','--');
ax3 = plot(nan,nan,'Color',cols.ptn_alm,'LineWidth',1.5,'LineStyle','-');
ax4 = plot(nan,nan,'Color',cols.ptn_m1,'LineWidth',1.5,'LineStyle','--');




%% fit parametric function to ACF

close all
clear tau

% ACF(x) = a * exp(-x/tau) * cos(x/(2*pi*t_period));
% where a is an overall amplitude
% tau is the decay timescale
% t_period is the oscillation period of the autocorrelation function
% x is time lag


% Define the function
% ACF = @(params, x) params(1) * exp(-x / params(2)) .* cos(x / (2 * pi * params(3)));
ACF = @(params, x) params(1) * ( exp(-x / params(2)) + params(3));

% data
x = linspace(lags(1),lags(end),numel(lags)); 
y = zscore(autocorr);  

% Initial guesses for parameters [a, tau, t_period]
x0 = [1 1 1];

% Define lower and upper bounds for parameters
lb = [0, 0, 0];  % All parameters should be positive
ub = [Inf, 20, Inf];  % No upper bounds
% lb = [-Inf -Inf -Inf];
% up = [Inf Inf Inf];



% % Fit the model
% f = figure;
% f.Position = [680   602   368   276];
% ax = prettifyAxis(gca);

% trust-region-reflective
options = optimset('Algorithm','trust-region-reflective','ScaleProblem','jacobian', 'FinDiffType','central','TolFun',1e-7);  % Optional: see iteration output
% options = optimset('Algorithm','levenberg-marquardt','ScaleProblem','jacobian', 'FinDiffType','central','FunctionTolerance',);  % Optional: see iteration output
% options = optimset('Algorithm','levenberg-marquardt', 'FinDiffType','central');  % Optional: see iteration output
for iunit = 1:size(autocorr,2)
    % if ~ptn_alm(iunit)
    %     continue
    % end

    % if ~ptn_m1(iunit) && ~ptn_alm(iunit)
    %     continue
    % elseif ptn_alm(iunit)
    %     c = cols.ptn_alm;
    % else
    %     c = cols.ptn_m1;
    % end
    
    best_params = lsqcurvefit(ACF, x0, x, y(:,iunit)',lb,ub, options);
    tau(iunit) = best_params(2);
    if tau(iunit)>30
        tau(iunit) = nan;
    end

    y_fit = ACF(best_params, x);

    % % Plot the fit
    % cla(ax)
    % hold on
    % plot(x.*params.dt*1000, y(:,iunit), 'b.', 'MarkerSize',20, 'DisplayName', 'Data'); hold on;
    % plot(x.*params.dt*1000, y_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit');
    % % plot(x.*params.dt*1000, y_fit, 'Color',c, 'LineWidth', 1, 'DisplayName', 'Fit');
    % % legend;
    % xlabel('lag (ms)');
    % ylabel('ACF(x)');
    % title(['Unit ' num2str(iunit)]);
    % grid on;
    % pause

end



%%

close all

ptn_alm = ismember(allquality,'tagged') & contains(allregion,'ALM');
ptn_m1 = ismember(allquality,'tagged') & contains(allregion,'M1');
un_alm = ~ismember(allquality,'tagged') & contains(allregion,'ALM');
un_m1 = ~ismember(allquality,'tagged') & contains(allregion,'M1');

cols = getColors;

cols.ptn_alm = [cols.ptn(1) cols.ptn(2)*1.3 cols.ptn(3)*1.3];
cols.ptn_m1 = cols.ptn./1.3;

cols.un_alm = cols.un.*1.3;
cols.un_m1 = cols.un./1.3;

% temp = tau.*params.dt*1000; % millisec
% temp = tau;
temp = interp1(lags,tau).*params.dt*10000;

xlims = [0 500];
% temp(temp<200 & (ptn_alm|un_alm)') = temp(temp<200 & (ptn_alm|un_alm)') + 30;
% temp(temp<200 & (ptn_m1|un_m1)') = temp(temp<200 & (ptn_m1|un_m1)') + 27;

normtype = 'pdf';

f = figure; 
f.Position = [502   558   663   291];
ax = prettifyAxis(subplot(1,2,1));
hold on;

h = histogram(temp(un_alm),'Normalization',normtype);
plot([nanmedian(temp(un_alm)) nanmedian(temp(un_alm))],[ax.YLim],'k--','Color',cols.un_alm,'LineWidth',2);
h.FaceColor = cols.un_alm; h.EdgeColor = 'none';
h = histogram(temp(ptn_alm),15,'Normalization',normtype);
plot([nanmedian(temp(ptn_alm)) nanmedian(temp(ptn_alm))],[ax.YLim],'k--','Color',cols.ptn_alm,'LineWidth',2);
h.FaceColor = cols.ptn_alm; h.EdgeColor = 'none';
xlim(xlims)
p_ranksum_alm = ranksum(temp(un_alm),temp(ptn_alm));
title(['ALM, p=' num2str(round(p_ranksum_alm,3))],'fontweight','normal','fontsize',12)
xlabel('Decay constant (ms)')
ylabel('pdf')
ax1 = plot(nan,nan,'Color','w');
ax2 = plot(nan,nan,'Color','w');
labels = {['median=',num2str(round(nanmedian(temp(un_alm)),2)) ' ms'],...
          ['PTN median=',num2str(round(nanmedian(temp(ptn_alm)),2)) ' ms']};
leg = legend([ax1 ax2],labels); leg.Box = 'off'; leg.FontSize=10;

ax = prettifyAxis(subplot(1,2,2));
hold on;
h = histogram(temp(un_m1),'Normalization',normtype);
plot([nanmedian(temp(un_m1)) nanmedian(temp(un_m1))],[ax.YLim],'k--','Color',cols.un_m1,'LineWidth',2);
h.FaceColor = cols.un_m1; h.EdgeColor = 'none';
h = histogram(temp(ptn_m1),15,'Normalization',normtype);
plot([nanmedian(temp(ptn_m1)) nanmedian(temp(ptn_m1))],[ax.YLim],'k--','Color',cols.ptn_m1,'LineWidth',2);
h.FaceColor = cols.ptn_m1; h.EdgeColor = 'none';
xlim(xlims)
p_ranksum_m1 = ranksum(temp(un_m1),temp(ptn_m1));
title(['tjM1, p=' num2str(round(p_ranksum_m1,3))],'fontweight','normal','fontsize',12)
ax1 = plot(nan,nan,'Color','w');
ax2 = plot(nan,nan,'Color','w');
labels = {['median=',num2str(round(nanmedian(temp(un_m1)),2)) ' ms'],...
          ['PTN median=',num2str(round(nanmedian(temp(ptn_m1)),2)) ' ms']};
leg = legend([ax1 ax2],labels); leg.Box = 'off'; leg.FontSize=10;



%% ranksum test - continuous distributions with equal medians?

p = ranksum(temp(un_alm),temp(ptn_alm))
p = ranksum(temp(un_m1),temp(ptn_m1))


%% permutation test

temp = interp1(lags,tau).*params.dt*10000;

nn = 20;
mu = nanmedian(temp(un_m1)) - nanmedian(temp(ptn_m1));
P = 0.5;
nboot = 100000;
for i = 1:nboot
    idx = randperm(nn*2);
    p1 = randsample(temp(un_m1),nn,true);
    p2 = randsample(temp(ptn_m1),nn,true);
    tt = cat(1,p1',p2');
    
    one = tt(idx(1:round(P*nn*2)));
    two = tt(idx(round(P*nn*2)+1:end));
    
    mu_(i) = nanmedian(one) - nanmedian(two);

end


pval = 1-sum(mu_<mu)/nboot




%% save
clear rez
rez.autocorr = autocorr;
rez.autocorr_shape = '(lags x neurons)';
rez.lags = lags';
rez.tau = (interp1(lags,tau).*params.dt*10000)';
rez.ptn_alm = ptn_alm;
rez.ptn_m1 = ptn_m1;
rez.un_alm = un_alm;
rez.un_m1 = un_m1;

save('results\autocorr\AutoCorrDecayConstant.mat','rez')



