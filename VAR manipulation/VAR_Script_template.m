% Make new clusgroup

data_masterdir = GetCurrentDataDir();
save_dir = data_masterdir;

% range_fish = GetFishRange();
range_fish = 1:18;


%% (Optional) Load fresh from file

% load(fullfile(save_dir,'VAR_new.mat'),'VAR');


%%
for i_fish = range_fish,

%     VAR = RenameClusGroup(VAR,i_fish,1,'Selection');

    VAR = RenameCluster(VAR,i_fish,[3,4],'Auto_OMR_M0.7');

%     VAR = CopyCluster(VAR,i_fish,[2,1],[6,1]);  
%     VAR = CopyCluster(VAR,i_fish,[2,2],[6,2]);  
%     VAR = CopyCluster(VAR,i_fish,[6,2],[2,1]);  
%     VAR = CopyCluster(VAR,i_fish,[6,1],[2,2]);  

%     VAR = MakeNewClusGroup(VAR,i_fish,'Other');

end

%%
save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');
% disp('saved updated ''VAR''');

%% save with backup

timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
arcmatfolder = fullfile(save_dir, 'arc mat');
save(fullfile(arcmatfolder,matname),'VAR','-v6');

save(fullfile(save_dir,'VAR_new.mat'),'VAR','-v6');


%% view all names
i_fish = 1;
nClusGroup = length(VAR(1).ClusGroupName);

VAR(i_fish).AllNames = [];
for i = 1:nClusGroup,
    head = VAR(i_fish).ClusGroupName(i);
    col = vertcat(head, {VAR(i_fish).ClusGroup{i}.name}');
    for j = 1:length(col),
        VAR(i_fish).AllNames{j,i} = col{j};
    end
end

