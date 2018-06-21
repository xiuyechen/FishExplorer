% Make new clusgroup

data_masterdir = GetCurrentDataDir();
save_dir = data_masterdir;

% range_fish = GetFishRange();
range_fish = 1:18;
global VAR;

%% (Optional) Load fresh from file

% load(fullfile(save_dir,'VAR_new.mat'),'VAR');


%%
range_fish = [1:18];
for i_fish = range_fish

    
%     VAR = RenameClusGroup(VAR,i_fish,6,'Auto_M0.7_woA');
%     VAR = DeleteCluster(VAR,i_fish,{12,1});
%     VAR = RenameCluster(VAR,i_fish,[4,1],'k20_hf_defS_CV1'); % 1 at a time

%     VAR = CopyCluster(VAR,i_fish,{6,1},{6,2});  

    VAR = MakeNewClusGroup(VAR,i_fish,'CV2_altframe');
% VAR = MoveCluster(VAR,i_fish,{10,2},{13,4});
%     VAR = MoveCluster(VAR,i_fish,{10,2:4},{13,1:3});


end
VAR = CompileVARnames(VAR);
% SaveVARwithBackup();
% save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');
% disp('saved updated ''VAR''');

%% save with backup

timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
arcmatfolder = fullfile(save_dir, 'arc mat');
save(fullfile(arcmatfolder,matname),'VAR','-v6');

save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');

