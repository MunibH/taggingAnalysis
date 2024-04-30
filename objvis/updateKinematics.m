function updateKinematics(~,~,fig)

h = guidata(fig);

cla(h.ax.kin)
ylim(h.ax.kin,"auto");

hold(h.ax.kin, 'on');

tm = h.obj.time;

trialid = h.trialid; % these conditions come from updatePSTH
cols = h.cols;
kin = h.kinstring;
ikin = h.kinix;

% heatmap of selected feature
alldat = [];
for i = 1:numel(trialid)
    this = squeeze(h.kin.dat(:,trialid{i},ikin));
    
    alldat = cat(2,alldat,this);

    if i>1
        nTrialsCond(i-1) = sum(cell2mat(cellfun(@numel,trialid(1:i-1),'uni',0 )));
    end
end

% alldat = abs(alldat);

mn = min(min(alldat));
mx = max(max(alldat));
imagesc(h.ax.kin,tm,1:size(alldat,2),alldat')
c = colorbar(h.ax.kin,"eastoutside","Limits",[mn,mx]);

% if ~isempty(h.kin_clim.String)
%     clim(h.ax.kin,eval(h.kin_clim.String))
% end


for i = 1:numel(nTrialsCond)
    plot(h.ax.kin,h.ax.kin.XLim,[nTrialsCond(i),nTrialsCond(i)],'--','LineWidth',1,'Color',[255, 89, 249]./255)
end



xlabel(h.ax.kin, ['Time from ' h.params.alignEvent ' (s)']);
ylabel(h.ax.kin, 'Trials');
title(h.ax.kin, ['Feat: ' strrep(kin,'_','-')]);

ylim(h.ax.kin,[h.ax.kin.YLim(1) size(alldat,2)+0.5])

guidata(fig, h);


end



