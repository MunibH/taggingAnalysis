function meta = allSessionMeta(meta,datapth)

% datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
% meta = [];

%% these sessions were before implementing tshift,tprime,bombcell

% H2 sessions
meta = loadJPV8(meta,datapth);
meta = loadJPV11(meta,datapth);
meta = loadJPV12(meta,datapth);
% NP2 sessions
meta = loadJPV13(meta,datapth);

%% these sessions were after implementing tshift,tprime,bombcell

end