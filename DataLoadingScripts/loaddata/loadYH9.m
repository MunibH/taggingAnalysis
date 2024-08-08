function meta = loadYH9(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'YH9';
meta(end).date = '2023-11-22';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = 'R_FN';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);


end
