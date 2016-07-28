function VAR = DeleteCluster(VAR,i_fish,sourceIDs)
% can delete several clusters at a time
nClusGroup = length(VAR(i_fish).ClusGroup);
clusgroup_list = intersect(1:nClusGroup,sourceIDs{1});

nClus = length(VAR(i_fish).ClusGroup{sourceIDs{1}});
clus_list = intersect(1:nClus,sourceIDs{2});

VAR(i_fish).ClusGroup{clusgroup_list}(clus_list) = [];
end
