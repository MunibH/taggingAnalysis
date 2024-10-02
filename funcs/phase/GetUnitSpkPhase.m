function spk = GetUnitSpkPhase(obj,cluid,lick)

%%

nLicks = lick.nLicks;

nProbes = numel(cluid);
nCluProbe = cell2mat(cellfun(@numel,cluid,'uni',0));

spk = cell(nProbes,1);
disp('Getting phase info for:')
for iprobe = 1:nProbes
    disp(['Probe ' num2str(iprobe)])
    for iunit = 1:nCluProbe(iprobe)
        if mod(iunit,10)==0
            disp(['Unit ' num2str(iunit) '/' num2str(nCluProbe(iprobe))])
        end

        thiscluid = cluid{iprobe}(iunit);
        clu = obj.clu{iprobe}(thiscluid);

        spk{iprobe}(iunit).cluid = thiscluid;

        spk{iprobe}(iunit).licktm = []; % time of spike after start of each lick
        spk{iprobe}(iunit).licknum = []; % lick number that spike occurred on, can use lick.trial to get trial number of that lick
        spk{iprobe}(iunit).phase = []; % phase of lick (0:2pi) that spike occurred


        for ilick = 1:nLicks
            
            % find which trial the lick occurred on
            on_trial_lick = lick.trial(ilick);
            % get frame times corresponding to the licks
            this_lick_trialtm = lick.frame_times{ilick};

            % get spikes within lick times on this trial
            on_trial_spike = clu.trial==on_trial_lick;
            in_time = (clu.trialtm>=this_lick_trialtm(1)) & (clu.trialtm<=this_lick_trialtm(end));

            useSpikes = (on_trial_spike & in_time);

            spks_during_lick = clu.trialtm(useSpikes);

            % align to lick onset
            spk{iprobe}(iunit).licktm = cat(1,spk{iprobe}(iunit).licktm , spks_during_lick - this_lick_trialtm(1));
            spk{iprobe}(iunit).licknum = cat(1,spk{iprobe}(iunit).licknum , ones(size(spks_during_lick))*ilick);

            % assign phase
            if ~isempty(spks_during_lick)
                thisphase = linspace(0,2*pi,numel(this_lick_trialtm));
                % find ix where spk times falls in frame_times, assign lick.phase as thisphase(ix)
                ix = findTimeIX(this_lick_trialtm,spks_during_lick);
                spk{iprobe}(iunit).phase = cat(1,spk{iprobe}(iunit).phase, thisphase(ix)');
            end

        end % ilick
    end % iunit
end  % iprobe



end