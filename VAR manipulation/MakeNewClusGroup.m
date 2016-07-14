function VAR = MakeNewClusGroup(VAR,i_fish,name)

nClusGroup = length(VAR(i_fish).ClusGroup);
i_ClusGroup = nClusGroup+1;

VAR(i_fish).ClusGroupName{i_ClusGroup} = name;   

% fill with holder:
temp = VAR(i_fish).ClusGroup{1}(1);
VAR(i_fish).ClusGroup{i_ClusGroup}(1).name = temp.name;
VAR(i_fish).ClusGroup{i_ClusGroup}(1).cIX_abs = temp.cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(1).gIX = temp.gIX;
VAR(i_fish).ClusGroup{i_ClusGroup}(1).numK = temp.numK;
    
end