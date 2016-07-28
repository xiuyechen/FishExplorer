function clusname = GetClusterName(VAR,i_fish,IDs)
clusname = VAR(i_fish).ClusGroup{IDs(1)}(IDs(2)).name;
end