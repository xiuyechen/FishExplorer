function [cIX,gIX,M_xyz_norm,M_xyz,cIX_abs,absIX] = GetDefaultClustersFromLoad(hfig,i_fish,option)

LoadFullFish(hfig,i_fish,-1);
absIX = getappdata(hfig,'absIX');
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
CellXYZ = getappdata(hfig,'CellXYZ');
if exist('option','var'),
    ClusterIDs = GetClusterIDs(option);
else
    ClusterIDs = GetClusterIDs();
end
[cIX,gIX] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2),absIX);
cIX_abs = absIX(cIX);
M_xyz_norm = CellXYZ_norm(cIX_abs,:);
M_xyz = CellXYZ(cIX_abs,:);