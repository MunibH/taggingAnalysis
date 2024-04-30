function [obj,params] = loadSessionData(meta,params, varargin)

% varargin only accepts one input right now -> params.behav_only. If 1,
% don't process clusters
if nargin > 2
    behav_only = varargin{1};
else
    behav_only = 0;
end

% load data obj
disp(["LOADING DATA OBJ: " meta.anm ' ' meta.date])
load(meta.datapth); % loads 'obj'

% delete tagging trials from obj if tagging session
if params.remove_tag
    obj = deleteTaggingTrials(obj);
end

% get event times
params.eventTimes = getEventTimes(obj.bp.ev,params.events,params.alignEvent);

% % add fake probe for testing
% meta.probe = [1 2];
% obj.clu{2} = obj.clu{1}(1:20);

% process data per specified probe
params.probe = meta.probe;
for prbix = 1:numel(params.probe)
    disp('______________________________________________________')
    disp(['Processing data for Probe ' num2str(params.probe(prbix))])
    disp(' ')

    prbnum = params.probe(prbix);

    [probeparams{prbix},probeobj{prbix}] = processData(obj,params,prbnum, behav_only);
end


if numel(params.probe) == 1
    params = probeparams{1};
    obj = probeobj{1};
else
    params = probeparams{1};
    params.cluid =  eval(eval(generateProbeString(numel(params.probe))));

    obj = probeobj{1};
    obj.psth = eval(eval(generatePSTHString(numel(params.probe))));
    obj.trialdat = eval(eval(generateTrialdatString(numel(params.probe))));
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')


% %%
% 
% % clean up sessparams and sessobj
% for sessix = 1:numel(meta)
%     params.trialid{sessix} = sessparams{sessix}.trialid;
% 
%     if numel(params.probe{sessix}) == 1
%         params.cluid{sessix} = sessparams{sessix,1}.cluid{params.probe{sessix}};
% 
%         objs(sessix) = sessobj{sessix,1};
% 
%         if behav_only
%             continue;
%         end
% 
%         objs(sessix).psth = objs(sessix).psth{params.probe{sessix}};
%         objs(sessix).trialdat = objs(sessix).trialdat{params.probe{sessix}};
%         % objs(sessix).presampleFR = objs(sessix).presampleFR{params.probe{sessix}};
%         % objs(sessix).presampleSigma = objs(sessix).presampleSigma{params.probe{sessix}};
%     elseif numel(params.probe{sessix}) == 2 % concatenate both probes worth of data
% 
%         params.cluid{sessix} = {sessparams{sessix,1}.cluid{params.probe{sessix}(1)}, sessparams{sessix,2}.cluid{params.probe{sessix}(2)} };
% 
%         objs(sessix) = sessobj{sessix,1};
% 
%         if behav_only
%             continue;
%         end
% 
%         objs(sessix).psth = cat(2, objs(sessix).psth{1}, sessobj{sessix,2}.psth{2});
%         objs(sessix).trialdat = cat(2, objs(sessix).trialdat{1}, sessobj{sessix,2}.trialdat{2});
%         % objs(sessix).presampleFR = cat(1, objs(sessix).presampleFR{1}, sessobj{sessix,2}.presampleFR{2});
%         % objs(sessix).presampleSigma = cat(1, objs(sessix).presampleSigma{1}, sessobj{sessix,2}.presampleSigma{2});
% 
%         % assign trialtm_aligned for 2nd probe
%         for i = 1:numel(objs(sessix).clu{2})
%             objs(sessix).clu{2}(i).trialtm_aligned = sessobj{sessix,2}.clu{2}(i).trialtm_aligned;
%         end
%     end
% end





% %% convert params to a struct array (just for convenience)
% 
% temp = params;
% clear params
% for sessix = 1:numel(meta)
%     temp2 = rmfield(temp,{'probe','trialid','cluid'});
%     temp2.probe = temp.probe{sessix};
%     temp2.trialid = temp.trialid{sessix};
%     temp2.cluid = temp.cluid{sessix};
%     params(sessix) = temp2;
% end
