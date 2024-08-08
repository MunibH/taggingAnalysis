function evtimes = getEventTimes(ev,events,alignev)

if strcmpi(alignev,'lastlick') || strcmpi(alignev,'firstlick')
    evtimes.(alignev) = 0;
    return
else
    align = mode(ev.(alignev));
end

for i = 1:numel(events)
    evtimes.(events{i}) = mode(ev.(events{i})) - align;
end




end