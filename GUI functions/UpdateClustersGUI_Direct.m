function UpdateClustersGUI_Direct(clusgroupID,clusID,i_fish)
% clusID = getappdata(hfig,'clusID');
% clusgroupID = getappdata(hfig,'clusgroupID');
% i_fish = getappdata(hfig,'i_fish');
global VAR;
ClusGroup = VAR(i_fish).ClusGroup{clusgroupID};

global hclusname hclusmenu;
set(hclusname,'String',ClusGroup(clusID).name);
menu = MakeNumberedMenu({ClusGroup.name});
set(hclusmenu,'String', menu,'Value',clusID+1);
end