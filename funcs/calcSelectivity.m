function sel = calcSelectivity(psth, cond2use)

sel = psth(:,:,cond2use(1)) - psth(:,:,cond2use(2));

end