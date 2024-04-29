function meta = loadJPV12(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV12';
meta(end).date = '2023-08-02';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'H2';
meta(end).region = 'L_ALM';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV12';
meta(end).date = '2023-08-05';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'H2';
meta(end).region = 'L_tjM1';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);


end
