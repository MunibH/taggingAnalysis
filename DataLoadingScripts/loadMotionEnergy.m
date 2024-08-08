function me = loadMotionEnergy(obj,meta,params,datapth)


% %------------------------------------------------------------
% % --first check to see if obj has an 'me' and has thresh--
% % -----------------------------------------------------------
objhasme = 0;
if isfield(obj,'me')
    me = obj.me;
    objhasme = 1;
end
% -----------------------------------------------------------
% -- if not, loadMotionEnergy, and save it it to obj.me--
% -----------------------------------------------------------

% find and load motion energy .mat file
if ~objhasme
    mepth = fullfile(datapth,'DataObjects',meta.anm);
    contents = dir(mepth);
    contents = contents(~[contents.isdir]);
    
    fns = {contents.name}';
    
    fn = patternMatchCellArray({contents.name}',{'motionEnergy',meta.date},'all');
    
    if numel(fn) == 1
        fn = fn{1};
    elseif numel(fn) == 2
        ix = ~contains(fn,'._'); % on windows, sometimes a file that starts with '._motionEnergy*' is found
        fn = fn{ix};
    else
        disp('UNABLE TO LOCATE A MOTION ENERGY FILE IN: ')
        disp(mepth)
        disp('      Continuing without motion energy as a feature')
        me.use = 0;
        me.data = nan;
        me.moveThresh = nan;
        return
    end
    
    temp = load(fullfile(mepth,fn));
    me = temp.me;
end

if iscell(me)
    temp = me;
    clear me
    me = struct();
    me.data = temp;
elseif isstruct(me.data)
    me.data = me.data.data;
end

% % % delete tagging trials (already doing this in deleteTaggingTrials
% % if params.remove_tag
% %     me.data = me.data(1:obj.bp.Ntrials); % first set of trials in motion energy are behavior trials
% % end



% -----------------------------------------------------------
% trim trial length (me.data contains motion energy for each time point in
% trial at 400 Hz). Want to align to params.alignEvent and want to put it
% in same dt as neural data
% -----------------------------------------------------------
taxis = obj.time + params.advance_movement;
alignTimes = obj.bp.ev.(params.alignEvent);
me.newdata = zeros(numel(obj.time),numel(me.data));
for trix = 1:obj.bp.Ntrials
    try
        me.newdata(:,trix) = interp1(obj.traj{1}(trix).frameTimes-0.5-alignTimes(trix),me.data{trix},taxis); % interp1(old_time,me,new_time);
    catch % if frameTimes doesn't exist or is full of NaNs - shouldn't be dummy data as we aren't using those sessions
        frameTimes = (1:size(obj.traj{1}(trix).ts,1)) ./ 400;
        me.newdata(:,trix) = interp1(frameTimes-0.5-alignTimes(trix),me.data{trix},taxis);
    end
end

% replace me.data with me.newdata
me.data = me.newdata;
me = rmfield(me,'newdata');
% fill nans with nearest value (there are some nans at the start of each
% trial)
me.data = fillmissing(me.data,'nearest');

% -------------------------------------------------------------------
% -- assign move time points as logical array same size as me.data --
% -------------------------------------------------------------------
% these values, except for YH come from 2-GMM fitting, see function below
% i've just rounded these to nearest whole number
if ~isfield(me,'moveThresh')
    if strcmp(meta.anm,'JPV8')
        me.moveThresh = 8;
    elseif strcmp(meta.anm,'JPV11') && strcmp(meta.date,'2023-06-16')
        me.moveThresh = 14;
    elseif strcmp(meta.anm,'JPV11') && strcmp(meta.date,'2023-06-21')
        me.moveThresh = 9;
    elseif strcmp(meta.anm,'JPV11') && strcmp(meta.date,'2023-06-22')
        me.moveThresh = 10;
    elseif strcmp(meta.anm,'JPV11') && strcmp(meta.date,'2023-06-23')
        me.moveThresh = 8;
    elseif strcmp(meta.anm,'JPV12') && strcmp(meta.date,'2023-08-02')
        me.moveThresh = 7;
    elseif strcmp(meta.anm,'JPV12') && strcmp(meta.date,'2023-08-05')
        me.moveThresh = 9;
    elseif strcmp(meta.anm,'JPV13') && strcmp(meta.date,'2023-10-03')
        me.moveThresh = 11;
    elseif strcmp(meta.anm,'JPV13') && strcmp(meta.date,'2023-10-04')
        me.moveThresh = 11;
    elseif strcmp(meta.anm,'JPV13') && strcmp(meta.date,'2023-10-05')
        me.moveThresh = 9;
    elseif strcmp(meta.anm,'MAH23') && strcmp(meta.date,'2024-07-03')
        me.moveThresh = 6;
    elseif strcmp(meta.anm,'MAH23') && strcmp(meta.date,'2024-07-05')
        me.moveThresh = 8;
    elseif strcmp(meta.anm,'MAH23') && strcmp(meta.date,'2024-07-06')
        me.moveThresh = 6;
    elseif strcmp(meta.anm,'MAH24') && strcmp(meta.date,'2024-06-11')
        me.moveThresh = 7;
    elseif strcmp(meta.anm,'MAH24') && strcmp(meta.date,'2024-06-12')
        me.moveThresh = 7;
    elseif strcmp(meta.anm,'MAH24') && strcmp(meta.date,'2024-06-14')
        me.moveThresh = 7;
    elseif strcmp(meta.anm,'MAH24') && strcmp(meta.date,'2024-06-15')
        me.moveThresh = 8;
    elseif strcmp(meta.anm,'YH11')
        me.moveThresh = 8;
    elseif strcmp(meta.anm,'YH9')
        me.moveThresh = 10.5;
    else
        error('New session, need to find movement threshold');
    end
end

% [~, ~, me.moveThresh] = detectMovementGMM(cell2mat(obj.me'));

me.move = me.data > (me.moveThresh);


% % baseline subtract 
% ps = [mode(obj.bp.ev.bitStart) mode(obj.bp.ev.sample-obj.bp.ev.bitStart)] - mode(obj.bp.ev.(params.alignEvent));
% for i = 1:numel(ps)
%     [~,psix(i)] = min(abs(obj.time - ps(i)));
% end
% psme = mean(me.data(psix(1):psix(2),:),1);
% me.data = me.data - psme;

end % loadMotionEnergy













