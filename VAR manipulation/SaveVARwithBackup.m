function SaveVARwithBackup(VAR)
save_dir = GetCurrentDataDir();
timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
arcmatfolder = fullfile(save_dir, 'arc mat');
save(fullfile(arcmatfolder,matname),'VAR','-v6');

save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');

end