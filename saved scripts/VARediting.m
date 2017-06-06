%% Load fish
global VAR;

range_fish = 1:18;

%%
for i_fish = range_fish
    %% load fish
%    LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    
    [cIX1,gIX1,~,absIX] = LoadCluster_Direct(i_fish,10,2);
    [cIX2,gIX2] = LoadCluster_Direct(i_fish,10,3);
    cIX = [cIX1;cIX2];
    gIX = [ones(size(cIX1));2*ones(size(cIX2))];
%     gIX = [gIX1;5-gIX2];
    
    name = 'HBO_LR';
    clusgroupID = 10;
    clusIDoverride = 5;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
    
end
%%     save VAR
SaveVARtoMat(hfig);