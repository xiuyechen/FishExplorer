function [cIX,gIX,numK] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX)
global VAR
ClusGroup = VAR(i_fish).ClusGroup{i_ClusGroup};
numK = ClusGroup(i_Cluster).numK;
gIX = ClusGroup(i_Cluster).gIX;
cIX_abs = ClusGroup(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
[~,cIX] = ismember(cIX_abs,absIX);
end