function plotEventTimes(ax,evtimes,varargin)

if nargin>2
    lc = varargin{1}; % line color
else
    lc = 'k';
end

if nargin>3
    setylim = varargin{2};
else
    setylim = true;
end

ylims = ax.YLim;

fnames = fieldnames(evtimes);
for i = 1:numel(fnames)
    this = evtimes.(fnames{i});
    plot(ax,[this,this],ylims,'--','Color',lc)
    % xline(ax,evtimes.(fnames{i}),'k--');
end
if setylim
    ylim(ax,ylims);
end

end