function [jawphase,jawamp] = getPhaseOfJaw(kin)

feats = {'jaw_ydisp_view1'}; % feature to get phase and amp of jaw from
featmask = ismember(kin.featLeg,feats);

jaw = kin.dat(:,:,featmask);

for trial = 1:size(jaw,2) % for each trial
    j = jaw(:,trial);
    [jawphase(:,trial),jawamp(:,trial)] = calculateJawPhase(j);
end


end
