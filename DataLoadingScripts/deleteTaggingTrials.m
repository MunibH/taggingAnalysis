% delete tagging trials from obj if tagging session
function obj = deleteTaggingTrials(obj)
if isfield(obj.bp, 'fns') % usually an optotagging session
    ix = find(cell2mat(cellfun(@(x) contains(x,'MasterProtocol'),obj.bp.fns,'uni',0)));
    behavMask = obj.bp.fidx==ix; % logical array - 1 if behav trial, 0 otherwise
    % obj.bp
    obj.bp.Ntrials = sum(behavMask);
    obj.bp.hit = obj.bp.hit(behavMask);
    obj.bp.miss = obj.bp.miss(behavMask);
    obj.bp.no = obj.bp.no(behavMask);
    obj.bp.early = obj.bp.early(behavMask);
    obj.bp.protocol.nums = obj.bp.protocol.nums(behavMask);
    obj.bp.bitRand = obj.bp.bitRand(behavMask);
    obj.bp.autowater = obj.bp.autowater(behavMask);
    obj.bp.R = obj.bp.R(behavMask);
    obj.bp.L = obj.bp.L(behavMask);
    if isfield(obj.bp,'autowaterBlock')
        obj.bp.autowaterBlock = obj.bp.autowaterBlock(behavMask);
    end
    obj.bp.autolearn = obj.bp.autolearn(behavMask);
    % obj.bp.ev
    obj.bp.ev.bitStart = obj.bp.ev.bitStart(behavMask);
    obj.bp.ev.sample = obj.bp.ev.sample(behavMask);
    obj.bp.ev.delay = obj.bp.ev.delay(behavMask);
    obj.bp.ev.goCue = obj.bp.ev.goCue(behavMask);
    obj.bp.ev.reward = obj.bp.ev.reward(behavMask);
    obj.bp.ev.lickL = obj.bp.ev.lickL(behavMask);
    obj.bp.ev.lickR = obj.bp.ev.lickR(behavMask);
    % obj.bp.stim
    obj.bp.stim.enable = obj.bp.stim.enable(behavMask);
    obj.bp.stim.num = obj.bp.stim.num(behavMask);
    % obj.bp.traj
    obj.traj{1} = obj.traj{1}(behavMask);
    obj.traj{2} = obj.traj{2}(behavMask);
    % obj.me
    obj.me = obj.me(1:obj.bp.Ntrials); % first set of trials is behav in motion energy cell array
    % obj.clu
    for iprb = 1:numel(obj.clu)
        for iclu = 1:numel(obj.clu{iprb})
            clu = obj.clu{iprb}(iclu);
            spkmask = ismember(clu.trial,find(behavMask));
            obj.clu{iprb}(iclu).trial = obj.clu{iprb}(iclu).trial(spkmask);
            % also have to subtract how many trials you removed from
            % beginning of session
            nTrials2Sub = find(behavMask==1,1,'first') - 1;
            obj.clu{iprb}(iclu).trial =  obj.clu{iprb}(iclu).trial - nTrials2Sub;
            obj.clu{iprb}(iclu).tm = obj.clu{iprb}(iclu).tm(spkmask);
            obj.clu{iprb}(iclu).trialtm = obj.clu{iprb}(iclu).trialtm(spkmask);
        end
    end
end
end