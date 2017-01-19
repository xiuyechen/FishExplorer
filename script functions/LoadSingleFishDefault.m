function [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange)

setappdata(hfig,'i_fish',i_fish);

%% load fish data
LoadFullFish(hfig,i_fish);

%% load the auto-clustering indices
if ~exist('ClusterIDs','var'),
    ClusterIDs = GetClusterIDs(); % Auto_M0.7_woA
end
i_ClusGroup = ClusterIDs(1);
i_Cluster = ClusterIDs(2); 

[cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);

setappdata(hfig,'cIX',cIX);
setappdata(hfig,'gIX',gIX);

if ~exist('stimrange','var'),
    [~,stimrange] = GetStimRange([],i_fish);
end
setappdata(hfig,'stimrange',stimrange); % need to UpdateTimeIndex if changed
UpdateTimeIndex(hfig);

M = getappdata(hfig,'M');
if nargout>3,
    stim = getappdata(hfig,'stim');
    behavior = getappdata(hfig,'behavior');
    M_0 = getappdata(hfig,'M_0');
end

end