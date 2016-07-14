function VAR = CopyCluster(VAR,i_fish,sourceIDs,destIDs)

Source = VAR(i_fish).ClusGroup{sourceIDs(1)}(sourceIDs(2));
VAR(i_fish).ClusGroup{destIDs(1)}(destIDs(2)) = Source;

end