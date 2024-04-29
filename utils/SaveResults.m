function SaveResults(pth,fn,varname,result)
% pth is where to save
% fn is file name


eval([varname, ' = result;']);

if ~exist(pth,'dir')
    mkdir(pth);
end

savepth = fullfile(pth,fn);

save(savepth,varname,'-v7.3')

disp(['Results saved ' varname ' to: ' savepth])

end