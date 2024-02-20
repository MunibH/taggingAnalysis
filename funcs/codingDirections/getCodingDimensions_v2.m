function rez = getCodingDimensions_v2(inp,obj,params)


cond2use = inp.cond2use;
cd_labels = inp.cd_labels;
cd_epochs = inp.cd_epochs;
cd_times = inp.cd_times;

% trials
for i = 1:numel(cond2use)
    trix{i} = params.trialid{cond2use(i)};
    nTrialsCond(i) = numel(trix{i});
end
% balance number of trials across conditions
nTrials = min(nTrialsCond);
for i = 1:numel(cond2use)
    trix{i} = randsample(trix{i},nTrials,false);
end

% get single trial data
data = cellfun(@(x) permute(obj.trialdat(:,:,x), [1 3 2]), trix, 'uni', 0); % (time,trials,units)

% standardize data
nTime = numel(obj.time);
nUnits = size(obj.trialdat,2);
data_reshaped = cellfun(@(x) reshape(x,nTime*nTrials,nUnits), data, 'uni', 0);
alldata = cat(1,data_reshaped{:});
mu = mean(alldata,1,'omitnan');
std_ = std(alldata,[],1,'omitnan');
for i = 1:numel(data)
    data_zscored{i} = (data_reshaped{i} - mu) ./ std_;
end
rez.trialdat = cellfun(@(x) reshape(x, nTime,nTrials,nUnits), data_zscored,'uni',0);


% compute psths
psth = cellfun(@(x) squeeze(mean(x,2)), rez.trialdat, 'uni', 0);
rez.psth = cat(3,psth{:});

% calc cds
align = mode(obj.bp.ev.(params.alignEvent));
rez.cd = zeros(size(rez.psth,2),numel(cd_labels)); % (neurons,numCDs)
for icd = 1:numel(cd_labels)
    % find time points to use
    edges = [cd_times{icd}(1) cd_times{icd}(2)];
    cdix = findTimeIX(obj.time,edges);
    cdix = cdix(1):cdix(2);
    % calculate coding direction
    rez.cd(:,icd) = calcCD(rez.psth,cdix);
end

% orthogonalize CDs
rez.cd_orth = gschmidt(rez.cd);

% calc projs
for i = 1:numel(rez.trialdat)
    rez.proj{i} = tensorprod(rez.trialdat{i},rez.cd_orth,3,1); % (time,trials,cds)
end


% calc VE





%%
% f = figure;
% ax = gca;
% hold on;
% for i = 1:nUnits
%     cla(ax)
%     plot(obj.time,mean(data{1}(:,:,i),2),'Color','b')
%     plot(obj.time,mean(data{2}(:,:,i),2),'Color','r')
%     pause
% end

% f = figure;
% ax = gca;
% hold on;
% for i = 1:nUnits
%     cla(ax)
%     plot(obj.time,psth{1}(:,i),'Color','b')
%     plot(obj.time,psth{2}(:,i),'Color','r')
%     pause
% end

%%


end % getCodingDimensions