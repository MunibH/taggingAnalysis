function meta = loadYH11(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'YH11';
meta(end).date = '2023-12-06';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = 'L_ALM';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'YH11';
meta(end).date = '2023-12-06';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 2;
meta(end).probeType = 'NP2';
meta(end).region = 'R_ALM';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);



end
