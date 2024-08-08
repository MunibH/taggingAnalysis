function meta = loadMAH24(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH24';
meta(end).date = '2024-06-11';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = [1,2];
meta(end).probeType = 'NP2';
meta(end).region = {'L_ALM','R_ALM'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH24';
meta(end).date = '2024-06-12';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = [1,2];
meta(end).probeType = 'NP2';
meta(end).region = {'L_ALM','R_ALM'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH24';
meta(end).date = '2024-06-14';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = {'R_tjM1'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH24';
meta(end).date = '2024-06-15';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = {'R_tjM1'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);



end
