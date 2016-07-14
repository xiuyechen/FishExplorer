function VAR = RenameCluster(VAR,i_fish,IDs,name)
VAR(i_fish).ClusGroup{IDs(1)}(IDs(2)).name = name;
end