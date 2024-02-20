function sel = calcPrefSelectivity_EachTimePoint(obj, params, cond2use, sm)


psth = obj.psth(:,:,cond2use);
% temp = permute(psth,[1 3 2]);
% sz = size(temp);
% temp = reshape(temp,sz(1)*sz(2),sz(3));
% temp = mySmooth(temp,sm,'zeropad');
% temp = zscore(temp);
% psth = reshape(temp,sz(1),sz(2),sz(3));
% psth = permute(psth,[1 3 2]);

sel = zeros(size(psth,1),size(psth,2));
for itime = 1:size(psth,1)
    dat = squeeze(psth(itime,:,:)); %(neurons,cond)
    [~,pref] = max(dat'); % preferred direction
    nonpref = 3 - pref; 
    for iclu = 1:size(psth,2)
        sel(itime,iclu) = dat(iclu,pref(iclu)) - dat(iclu,nonpref(iclu));
    end
end

end