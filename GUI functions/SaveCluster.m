function SaveCluster(hfig,cIX,gIX,name,clusgroupID)
absIX = getappdata(hfig,'absIX');
i_fish = getappdata(hfig,'i_fish');

global VAR;
ClusGroup = VAR(i_fish).ClusGroup{clusgroupID};

clusID = numel(ClusGroup)+1;
ClusGroup(clusID).name = name;

ClusGroup(clusID).cIX_abs = absIX(cIX);
ClusGroup(clusID).gIX = gIX;

U = unique(gIX);
numU = length(U);
ClusGroup(clusID).numK = numU;

VAR(i_fish).ClusGroup{clusgroupID} = ClusGroup;

UpdateClustersViewGUI(hfig);
disp('cluster saved');
end