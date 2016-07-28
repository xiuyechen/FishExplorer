function VAR = CopyCluster(VAR,i_fish,sourceIDs,destIDs)
nClus = length(VAR(i_fish).ClusGroup{sourceIDs{1}});
clus_list = intersect(1:nClus,sourceIDs{2});

Source = VAR(i_fish).ClusGroup{sourceIDs{1}}(clus_list);

dest_clus_list = destIDs{2}(1:length(clus_list));
VAR(i_fish).ClusGroup{destIDs{1}}(dest_clus_list) = Source;

end