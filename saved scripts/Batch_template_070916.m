% batch input data setup:
isFullData = 1;
data_masterdir = GetCurrentDataDir();

const_ClusGroup = 2;
const_Cluster = 1; % This is all cells
M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish =  1:18; % range_fish = GetFishRange();

%% custom params here:


%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    fishset = M_fishset(i_fish);

    LoadFullFish(hfig,i_fish,isFullData);
    absIX = getappdata(hfig,'absIX');
    
    %% Cluster indexing
    i_ClusGroup = const_ClusGroup;% M_ClusGroup(i);
    i_Cluster = const_Cluster;% M_Cluster(i);
    cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);

    %% get time index
    timelists = getappdata(hfig,'timelists');
%     timelists_names = getappdata(hfig,'timelists_names'); % use to screen
    stimrange = M_stimrange(i);
    
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);

    %% get data matrix
    M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    %M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
        
    %% ------custom code here---------

%         gIX = Kmeans_Direct(M,numK);
        

    %% (optional) save cluster
%     name = 'test';
%     clusgroupID = 2;
%     clusIDoverride = ; %Saved in Group2, Clusters 3&4
%     SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);

end