function ax = prettifyAxis(ax,varargin)

tl = 2; % tick length multiplier
fs = 13; % axis fontsize
lw = 2; % axis linewidth

if nargin > 1
    keys = varargin(1:2:end);
    values = varargin(2:2:end);

    if ~iscell(keys)
        k = keys;
        v = values;
        eval([k '= v;']);
    else
        for i = 1:numel(keys)
            k = keys{i};
            v = values{i};
            eval([k '= v;']);
        end
    end
end

% make axes black
set(groot, 'DefaultAxesXColor', [0,0,0], ...
'DefaultAxesYColor', [0,0,0], ...
'DefaultAxesZColor', [0,0,0]);

% change line thicknesses
ax.LineWidth = lw;

% change tick direction to outside
ax.TickDir = 'out';

% change tick size
ax.TickLength = ax.TickLength .* tl;

% change axis font size
ax.FontSize = fs;


end