function meta = loadJPV8(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV8';
meta(end).date = '2023-04-06';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'H2';
meta(end).region = {'L_ALM'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);


end
