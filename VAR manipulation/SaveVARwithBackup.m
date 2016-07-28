function SaveVARwithBackup()
global VAR; 
save_dir = GetCurrentDataDir();
timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
arcmatfolder = fullfile(save_dir, 'arc mat');
save(fullfile(arcmatfolder,matname),'VAR','-v6');

save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');

VAR = CompileVARnames(VAR);

disp('saved VAR with backup');
end