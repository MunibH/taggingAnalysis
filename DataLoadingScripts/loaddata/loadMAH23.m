function meta = loadMAH23(meta,datapth)

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH23';
meta(end).date = '2024-07-03';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = {'L_tjM1'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

% meta(end+1).datapth = datapth;
% meta(end).anm = 'MAH23';
% meta(end).date = '2024-07-04'; % only 1 tagged unit but it failed collision test - not using this session
% meta(end).datafn = findDataFn(meta(end));
% meta(end).probe = 1;
% meta(end).probeType = 'NP2';
% meta(end).region = 'L_tjM1';
% meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH23';
meta(end).date = '2024-07-05';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = {'R_tjM1'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);

meta(end+1).datapth = datapth;
meta(end).anm = 'MAH23';
meta(end).date = '2024-07-06';
meta(end).datafn = findDataFn(meta(end));
meta(end).probe = 1;
meta(end).probeType = 'NP2';
meta(end).region = {'R_ALM'};
meta(end).datapth = fullfile(meta(end).datapth,'DataObjects',meta(end).anm,meta(end).datafn);


end
