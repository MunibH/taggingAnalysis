function meta = allSessionMeta(meta,datapth)

% datapth = 'C:\Users\munib\Documents\Economo-Lab\data';
% meta = [];

meta = loadJPV8(meta,datapth);
meta = loadJPV11(meta,datapth);
meta = loadJPV12(meta,datapth);
meta = loadJPV13(meta,datapth);

end