function kin = getKinematics(obj,me,params)

% ----------------------------------------------
% -- Check if obj already has kin --
% ----------------------------------------------
if isfield(obj,'kin')
    if isfield(params,'reload_kin')
        if params.reload_kin
            % get kinematics
        else
            return % kin = obj.kin
        end
    end
    return % kin = obj.kin
else
    % get kinematics
end


% ----------------------------------------------
% -- Kinematics --
% (xdisp,ydisp,xvel,yvel) for each dlc feature
% ----------------------------------------------
[kin.dat,kin.featLeg,kin.nans] = getKinematicsFromVideo(obj,params);

% ----------------------------------------------
% -- Tongue angle and length --
% ----------------------------------------------
if ismember('tongue',params.traj_features{1})
    [ang,len] = getLickAngleAndLength(kin.featLeg,kin.dat,kin.nans);
    kin.featLeg{end+1} = 'tongue_angle';
    kin.featLeg{end+1} = 'tongue_length';

    kin.dat(:,:,end+1) = ang;
    kin.dat(:,:,end+1) = len;
end


% ----------------------------------------------
% -- Phase of jaw --
% ----------------------------------------------
[jawphase,jawamp] = getPhaseOfJaw(kin); %calculateJawPhase(jawpos(:,ct));
kin.featLeg{end+1} = 'jaw_phase';
kin.featLeg{end+1} = 'jaw_amp';

kin.dat(:,:,end+1) = jawphase;
kin.dat(:,:,end+1) = jawamp;

% ----------------------------------------------
% -- Motion energy --
% ----------------------------------------------
kin.featLeg{end+1} = 'motion_energy';

kin.dat(:,:,end+1) = me.data;

% % ----------------------------------------------
% % -- Standardize kinematics --
% % ----------------------------------------------
% for featix = 1:size(kin.dat,3)
%     temp = kin.dat(:,:,featix);
%     temp2 = temp(:);
%     kin.dat_std(:,:,featix) = (temp - mean(temp2,'omitnan')) ./ std(temp2,'omitnan'); 
% end

% % ----------------------------------------------
% % -- DIMENSIONALITY REDUCTION --
% % ----------------------------------------------
% % many of the video features will be highly correlated, so we will perform PCA/FA
% % on the matrix of features to reduce the dimensionality to a set of factors that
% % best explain the movement captured by the video recordings
% kin.dat_reduced = reduceDimensionVideoFeatures(kin.dat,params.feat_varToExplain);

kin = rmfield(kin,'nans');

end

