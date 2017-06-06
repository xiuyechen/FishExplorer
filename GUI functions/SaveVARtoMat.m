function SaveVARtoMat(hfig)
save_dir = getappdata(hfig,'save_dir');

% copy of VAR files will be saved into this subfolder:
arcmatfolder = fullfile(save_dir, 'arc mat');
if ~exist(arcmatfolder, 'dir')
    mkdir(arcmatfolder);
end

global VAR; %#ok<NUSED>
timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
save(fullfile(arcmatfolder,matname),'VAR','-v6');

% also save the current VAR file
save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');
disp('saved both to workspace and .mat');

end