% batch input data setup:
isFullData = 1;
data_masterdir = GetCurrentDataDir();

const_ClusGroup = 2;
const_Cluster = 1; % This is all cells
M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish =  1:18; % range_fish = GetFishRange();

%% custom params here:
numK = 20;

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
    % Load fish
    LoadFullFish(hfig,i_fish,isFullData);

    %% 1. 
    % setup
    absIX = getappdata(hfig,'absIX');
    fishset = M_fishset(i_fish);
    i_ClusGroup = const_ClusGroup;% M_ClusGroup(i);
    i_Cluster = const_Cluster;% M_Cluster(i);
    stimrange = M_stimrange{i_fish};
    timelists = getappdata(hfig,'timelists');
    
    % Load cluster data
    cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    %M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
        
    % ------custom code here---------      
    gIX = Kmeans_Direct(M,numK);
    
    % save cluster
    name = 'k20_full_defS';
    clusgroupID = 3;
    clusIDoverride = 1;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);

	%% 2. 
    % setup
    i_ClusGroup = 2;
    i_Cluster = 2;
    
    % Load cluster data
    cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
        
    % ------custom code here---------      
    gIX = Kmeans_Direct(M,numK);
    
    % save cluster
    name = 'k20_half_defS';
    clusgroupID = 3;
    clusIDoverride = 2;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
   	%% 3. 
    if ismember(i_fish,[8:15,17:18]),
        % setup
        i_ClusGroup = 2;
        i_Cluster = 1;
        stimrange = 1;
        
        % Load cluster data
        cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
        tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
        M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
        
        % ------custom code here---------
        gIX = Kmeans_Direct(M,numK);
        
        % save cluster
        name = 'k20_full_PT';
        clusgroupID = 3;
        clusIDoverride = 3;
        SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
        
        %%
        % setup        
        i_ClusGroup = 2;
        i_Cluster = 1;
        stimrange = 2;
        
        % Load cluster data
        cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
        tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
        M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
        
        % ------custom code here---------
        gIX = Kmeans_Direct(M,numK);
        
        % save cluster
        name = 'k20_full_OMR';
        clusgroupID = 3;
        clusIDoverride = 4;
        SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    end
end

SaveVARwithBackup();
