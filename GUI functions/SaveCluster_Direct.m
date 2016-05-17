function clusID = SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride)
% absIX = getappdata(hfig,'absIX');
% i_fish = getappdata(hfig,'i_fish');

global VAR;
ClusGroup = VAR(i_fish).ClusGroup{clusgroupID};

if exist('clusIDoverride','var'),
    clusID = clusIDoverride;
else
    clusID = numel(ClusGroup)+1;
end
ClusGroup(clusID).name = name;

ClusGroup(clusID).cIX_abs = absIX(cIX);
ClusGroup(clusID).gIX = gIX;

U = unique(gIX);
numU = length(U);
ClusGroup(clusID).numK = numU;

VAR(i_fish).ClusGroup{clusgroupID} = ClusGroup;

% UpdateClustersViewGUI(hfig);
disp('cluster saved');
end