% Multi-fish info

% save all CellXYZ in .mat
range_fish = 1:18;
emptycells = cell(length(range_fish),1);
CellXYZ_multi = emptycells;
CellXYZ_norm_multi = emptycells;
cIX_Auto07_multi = emptycells;
gIX_Auto07_multi = emptycells;
cIX_Auto05_multi = emptycells;
gIX_Auto05_multi = emptycells;
Clusmean_Auto07_multi = emptycells;
Clusmean_Auto05_multi = emptycells;

% Centroids_Auto07_multi = emptycells;
% Centroids_Auto05_multi = emptycells;
Behavior_multi = emptycells;
Stim_multi = emptycells;
numClus_multi = zeros(length(range_fish),1);
%
for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    % i_fish = 8;
    disp(i_fish);
    
    isFullData = true;
    LoadFullFish(hfig,i_fish,isFullData);
    CellXYZ_multi{i_fishnum} = getappdata(hfig,'CellXYZ');
    CellXYZ_norm_multi{i_fishnum} = getappdata(hfig,'CellXYZ_norm');
    
    %% load the auto-clustering results (indices only, not stim-range confusion)
    i_ClusGroup = 6;
    i_Cluster = 1;
    
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
    cIX_Auto07_multi{i_fishnum} = cIX;
    gIX_Auto07_multi{i_fishnum} = gIX;
    numClus_multi(i_fishnum) = max(gIX);
     
    %% load cluster time-series data
    M_stimrange = GetStimRange();    
    stimrange = M_stimrange{i_fish};    
    timelists = getappdata(hfig,'timelists');
    fishset = getappdata(hfig,'fishset');

    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    [M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    
    Clusmean_Auto07_multi{i_fishnum} = FindClustermeans(gIX,M);
    Behavior_multi{i_fishnum} = behavior;
    Stim_multi{i_fishnum} = stim;
    
    %% (another set)
    i_ClusGroup = 7;
    i_Cluster = 1;
    
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
    cIX_Auto05_multi{i_fishnum} = cIX;
    gIX_Auto05_multi{i_fishnum} = gIX;
    
    M_stimrange = GetStimRange();    
    stimrange = M_stimrange{i_fish};    
    fishset = getappdata(hfig,'fishset');

    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    
    Clusmean_Auto05_multi{i_fishnum} = FindClustermeans(gIX,M);
end
currdir = GetCurrentDataDir();
save(fullfile(currdir,'BasicMultiFishInfo.mat'),'CellXYZ_multi','CellXYZ_norm_multi',...
    'cIX_Auto07_multi','gIX_Auto07_multi','cIX_Auto05_multi','gIX_Auto05_multi',...
    'Clusmean_Auto07_multi','Clusmean_Auto05_multi');