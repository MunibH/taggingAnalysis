function ax = prettifyAxis(ax,varargin)

if nargin > 1
    tl = varargin{1}; % tick length multiplier
else 
    tl = 2;
end

if nargin > 2
    fs = varargin{2}; % axis font size
else 
    fs = 13;
end


% make axes black
set(groot, 'DefaultAxesXColor', [0,0,0], ...
'DefaultAxesYColor', [0,0,0], ...
'DefaultAxesZColor', [0,0,0]);

% change line thicknesses
ax.LineWidth = 1;

% change tick direction to outside
ax.TickDir = 'out';

% change tick size
ax.TickLength = ax.TickLength .* tl;

% change axis font size
ax.FontSize = fs;


end