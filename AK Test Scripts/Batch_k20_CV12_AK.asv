% batch
isFullData = 1;
numK = 20;
data_masterdir = GetCurrentDataDir();

range_fish =  2:2; 
% M_ClusGroup = [2,2,2,2];
% M_Cluster = [1,1,1,1];
const_ClusGroup = 2;
const_Cluster = 2; % This is all cells
% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]; 
M_stim = [1]; % PT only, for now
%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);

    LoadFullFish(hfig,i_fish,isFullData);
    absIX = getappdata(hfig,'absIX');
    
    %% Cluster indexing
    i_ClusGroup = const_ClusGroup;% M_ClusGroup(i);
    i_Cluster = const_Cluster;% M_Cluster(i);
    cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);

    %% partitions for CV
    timelists = getappdata(hfig,'timelists');
    timelists_names = getappdata(hfig,'timelists_names');
    periods = getappdata(hfig,'periods');
    %if length(periods)>1,
        timelistsCV = cell(length(M_stim),2);
        k_stim = 1;
         for k_stim = 1:length(M_stim), % :3
            i_stim = M_stim(k_stim);
            TL = timelists{i_stim};
            period = periods(i_stim);
            nrep = size(TL,2)/periods(i_stim); % integer
            n = floor(nrep/2);
            timelistsCV{k_stim,1} = TL(1):TL(n*period);
            timelistsCV{k_stim,2} = TL(1+n*period):TL(2*n*period);
         end
   % end
    
    %         for k_stim = 1:length(M_stim), % :3
 %%
    for k = 1:2,% CV halves
        tIX = timelistsCV{k_stim,k};
        M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX); 
        %M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
        

        gIX = Kmeans_Direct(M,numK);
        
        % save        
        name = ['k20_' num2str(k)];
        clusgroupID = 2;
        clusIDoverride = k+2; %Saved in Group2, Clusters 3&4
        SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
        %     end
    end
end