function [jawphase,jawamp] = getPhaseOfJaw(obj,kin,dt)

% feats = {'jaw_ydisp_view1'}; % feature to get phase and amp of jaw from
% featmask = ismember(kin.featLeg,feats);
% 
% jaw = kin.dat(:,:,featmask);

traj = obj.traj{1}; % side view
feat = 'jaw';
ifeat = ismember(traj(1).featNames,feat);

%%
% f = figure;
% ax = gca;
% hold on;
for trial = 1:numel(traj) % for each trial
    % yyaxis right
    % cla(ax)
    % yyaxis left
    % cla(ax)
    % % % % j = jaw(:,trial);
    jaw = traj(trial).ts(:,2,ifeat); % jaw, y coord
    [jawphase{trial},jawamp{trial},analytic{trial},jaw_filt{trial}] = calculateJawPhase(jaw,400);
    % yyaxis left
    % plot(jaw,'-','LineWidth',1.3)
    % ylabel('jaw pos')
    % yyaxis right
    % plot(jawphase{trial},'-','LineWidth',1.4)
    % ylabel('jaw phase (rad)')
    % pause
end
%%

end
