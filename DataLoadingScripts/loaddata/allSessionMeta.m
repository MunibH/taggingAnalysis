function meta = allSessionMeta(meta,datapth)

% H2 sessions
meta = loadJPV8(meta,datapth);  % 1 session
meta = loadJPV11(meta,datapth); % 4 sessions
meta = loadJPV12(meta,datapth); % 2 sessions
meta = loadJPV13(meta,datapth); % 3 sessions
% NP2 sessions
meta = loadMAH23(meta,datapth); % 3 sessions
meta = loadMAH24(meta,datapth); % 4 sessions (2 dual-probe)

end