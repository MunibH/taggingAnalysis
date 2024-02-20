function meta = loadJPV11(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV11';
meta(end).date = '2023-06-16';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).region = 'L_tjM1';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV11';
meta(end).date = '2023-06-21';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).region = 'R_tjM1';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV11';
meta(end).date = '2023-06-22';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).region = 'R_ALM';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'JPV11';
meta(end).date = '2023-06-23';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).region = 'R_tjM1';
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

end
