function mysavefig(f,pth,fn, varargin)
% f is fig handle
% pth is where to save
% fn is file name

if nargin > 3
    saveSVG = varargin{1};
else
    saveSVG = 0;
end

if ~exist(pth,'dir')
    mkdir(pth)
end

savepth = fullfile(pth,fn);

savefig(f,savepth)
saveas(f,savepth,'png')
if saveSVG
    saveas(f,savepth,'svg')
end
%     saveas(f,fullfile(pth,fn),'epsc')

disp(['Figure saved to: ' savepth])

end












