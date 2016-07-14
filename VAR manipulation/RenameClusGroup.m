function VAR = RenameClusGroup(VAR,i_fish,i_ClusGroup,name)
VAR(i_fish).ClusGroupName{i_ClusGroup} = name;
end