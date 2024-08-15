function [seq,time] = getSeqLickTriggered(phase,licknum,dt,sm)

% edges = params.tmin:params.dt:params.tmax;
% obj.time = edges + params.dt/2;
% obj.time = obj.time(1:end-1);
edges = 0:dt:(2*pi);
time = edges + dt/2;
time = time(1:end-1);

nLicks = max(licknum);

% get single trial data
seq = zeros(numel(time),nLicks);

for ilick = 1:nLicks

    spkix = ismember(licknum, ilick);

    % if no spikes found for current trial (j), move on to
    % next trial
    if all(~spkix)
        continue
    end

    N = histc(phase(spkix), edges);
    N = N(1:end-1);
    if size(N,2) > size(N,1)
        N = N'; % make sure N is a column vector
    end

    seq(:,ilick) = mySmooth(N./dt,sm,'reflect');

end


end