function VAR = CompileVARnames(VAR)
%% view all names
range_fish = 1:length(VAR);
for i_fish = range_fish;
    nClusGroup = length(VAR(i_fish).ClusGroupName);
    
    VAR(i_fish).AllNames = [];
    for i = 1:nClusGroup,
        head = VAR(i_fish).ClusGroupName(i);        
        col = vertcat(head, {'-----'} ,{VAR(i_fish).ClusGroup{i}.name}');
        for j = 1:length(col),
            if j>2,
                VAR(i_fish).AllNames{j,i} = [num2str(j-2),': ',col{j}];
            else                
                VAR(i_fish).AllNames{j,i} = col{j};
            end
        end
    end
end
end