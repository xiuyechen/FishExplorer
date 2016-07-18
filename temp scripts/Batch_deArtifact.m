% batch input data setup:
isFullData = 1;
data_masterdir = GetCurrentDataDir();

const_ClusGroup = 2;
const_Cluster = 1; % This is all cells
M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish =  1:18; % range_fish = GetFishRange();

%% custom params here:
% numK = 20;

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
    % Load fish
    LoadFullFish(hfig,i_fish,isFullData);
    
    absIX = getappdata(hfig,'absIX');
    fishset = M_fishset(i_fish);
    i_ClusGroup = 3;
    i_Cluster = 1;
%     stimrange = M_stimrange{i_fish};
%     timelists = getappdata(hfig,'timelists');
%     tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%     M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    %M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
    
    % Load cluster data
    [cIX_last,gIX_last] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
    
%     cIX_last = getappdata(hfig,'cIX');
%     gIX_last = getappdata(hfig,'gIX');
    CellXYZ = getappdata(hfig,'CellXYZ');
%     absIX = getappdata(hfig,'absIX');
    
    [cIX,~] = ArtifactAnalysis(cIX_last,gIX_last,CellXYZ,absIX);
    if isempty(cIX),
        disp('No artifacts found');
    end
    
    % set-difference
    [cIX,ia] = setdiff(cIX_last,cIX);
    gIX = gIX_last(ia);
    numK = length(unique(gIX));
    
    [gIX, numK] = SqueezeGroupIX(gIX);
    
        % save cluster
    name = 'Auto_defS_woA_Master0.5';
    clusgroupID = 7;
    clusIDoverride = 1;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
end