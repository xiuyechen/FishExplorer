function UpdateClustersViewGUI(hfig)
clusID = getappdata(hfig,'clusID');
clusgroupID_view = getappdata(hfig,'clusgroupID_view');
i_fish = getappdata(hfig,'i_fish');
global VAR;
ClusGroup = VAR(i_fish).ClusGroup{clusgroupID_view};

global hclusname hclusmenu;
set(hclusname,'String',ClusGroup(clusID).name);
menu = MakeNumberedMenu({ClusGroup.name});
set(hclusmenu,'String', menu,'Value',clusID+1);
end