% batch input data setup:
tic
isFullData = 1;
data_masterdir = GetCurrentDataDir();

const_ClusGroup = 2;
const_Cluster = 2; % This is 50%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cells
M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish =  8:18; % range_fish = GetFishRange();

%% custom params here:
numK1 = 20;
isWkmeans = false;
isMakeFoxels = false;
masterthres = 0.5;
clusParams = struct('merge',masterthres,'cap',masterthres,'reg1',masterthres,...
    'reg2',masterthres,'minSize',10,'k1',numK1);

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
    i_count = 1;
    i_ClusGroup = 2;
    i_Cluster = 2+i_count;
    
    stimrange = M_stimrange{i_fish};
    timelists = getappdata(hfig,'timelists');
    
    % Load cluster data
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    %     M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
    
    % ------custom code here---------
    [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
    
    % save cluster
    name = 'Foxels_defStim_reg0.7';
    clusgroupID = 5;
    clusIDoverride = i_count;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
    % ------custom code here---------
    [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
    
    % save cluster
    name = 'Auto_defStim_Master0.7';
    clusgroupID = 3;
    clusIDoverride = i_count;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
    %% 2.
    % setup
    i_count = 2;
    i_ClusGroup = 2;
    i_Cluster = 2+i_count;
    
    % Load cluster data
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
    
    % ------custom code here---------
    [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
    
    % save cluster
    name = 'Foxels_half_defStim_reg0.7';
    clusgroupID = 5;
    clusIDoverride = i_count;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
    % ------custom code here---------
    [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
    
    % save cluster
    name = 'Auto_half_defStim_Master0.7';
    clusgroupID = 3;
    clusIDoverride = i_count;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
end

%% 3.
% for i_fish = 8:18,
%     disp(i_fish);
%     
%     % Load fish
%     LoadFullFish(hfig,i_fish,isFullData);
% 
%     % setup
%     absIX = getappdata(hfig,'absIX');
%     fishset = M_fishset(i_fish);
%     timelists = getappdata(hfig,'timelists');
% 
%     if ismember(i_fish,[8:15,17:18]),
%         % setup
%         i_count = 3;
%         i_ClusGroup = 2;
%         i_Cluster = 2+i_count;
%         stimrange = 1;
%         
%         % Load cluster data
%         [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%         tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%         M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
%         
%         % ------custom code here---------
%         [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
%         
%         % save cluster
%         name = 'Foxels_PT_defStim_reg0.7';
%         clusgroupID = 5;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%         
%         % ------custom code here---------
%         [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
%         
%         % save cluster
%         name = 'Auto_PT_defStim_Master0.7';
%         clusgroupID = 3;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%         
%         
%         %%
%         % setup
%         i_count = 4;
%         i_ClusGroup = 2;
%         i_Cluster = 2+i_count;
%         stimrange = 2;
%         
%         % Load cluster data
%         [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%         tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%         M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
%         
%         % ------custom code here---------
%         [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
%         
%         % save cluster
%         name = 'Foxels_OMR_defStim_reg0.7';
%         clusgroupID = 5;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%         
%         % ------custom code here---------
%         [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
%         
%         % save cluster
%         name = 'Auto_OMR_defStim_Master0.7';
%         clusgroupID = 3;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%     end
%     
% end
toc

