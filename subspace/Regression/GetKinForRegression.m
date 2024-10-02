function kin = GetKinForRegression(obj,me,par)
feats2use = {'tongue_xdisp_view1','tongue_ydisp_view1','jaw_xdisp_view1','jaw_ydisp_view1', ...
    'top_tongue_xdisp_view2','top_tongue_ydisp_view2','motion_energy'};

kin = getKinematics(obj, me, par);
temp = kin.dat(:,:,ismember(kin.featLeg,feats2use));

for i = 1:numel(feats2use)
    ts = temp(:,:,i);
    if contains(feats2use{i},'tongue')
        mu = mean(mode(ts));
    else
        mu = nanmean(ts,[1 2]);
    end
    sigma = nanstd(ts,0,[1,2]);

    ts = (ts - mu) ./ sigma;
    kin_zscored(:,:,i) = ts;
end

kin = kin_zscored;

end