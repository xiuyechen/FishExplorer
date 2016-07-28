% batch input data setup:
isFullData = -1;
data_masterdir = GetCurrentDataDir();

M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish =  1:18; % range_fish = GetFishRange();

global VAR;
%% custom params here:
% numK = 20;

%%
range_fish =  8:18;

%     if i_Cluster>1,
%         range_fish =  8:18;
%     end
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    % Load fish
    LoadFullFish(hfig,i_fish,isFullData);
    
    absIX = getappdata(hfig,'absIX');
    CellXYZ = getappdata(hfig,'CellXYZ');
    
    for i_Cluster = 2:3,
        i_ClusGroup = 8;
        [cIX_last,gIX_last] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
        
        clusname = GetClusterName(VAR,i_fish,[i_ClusGroup,i_Cluster]);

        % custom code:
        [cIX,~] = ArtifactAnalysis(cIX_last,gIX_last,CellXYZ,absIX);
        if isempty(cIX),
            disp('No artifacts found');
        end
        
        % set-difference
        [cIX,ia] = setdiff(cIX_last,cIX);
        gIX = gIX_last(ia);
        
        [gIX, numK] = SqueezeGroupIX(gIX);
        
        % save cluster
        name = [clusname,'_woA'];
        clusgroupID = 6;
        clusIDoverride = i_Cluster;
        SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    end
end
SaveVARwithBackup();