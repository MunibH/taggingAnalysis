function cd = calcCD(psth,times)
% calculate a coding direction given:
% psth - (time,neurons,cond2use) (cond2use = 2 dimensions)
% times - indices or logical array of size(psth,1)=time
%           specifies time points to use when calculating CD

    mu = squeeze(mean(psth(times,:,:),1));

    sd = squeeze(std(psth(times,:,:),[],1));

    cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
    cd(isnan(cd)) = 0;
%     cd = cd ./ norm(cd);
    cd = cd./sum(abs(cd)); % (ncells,1)
end