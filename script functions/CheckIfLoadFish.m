function [M_fishrange_im,fishrange_load] = CheckIfLoadFish(M_fishrange,M_ClusterIDs,M_stimrange,M_fishrange_im)
if ~exist('M_fishrange_im','var')
    M_fishrange_im = M_fishrange;
end

global VAR;
for i_set = 1:length(M_fishrange)
    for i_fish = M_fishrange{i_set}
        
        % check M_fishrange
        flag1 = ismember(i_fish,cell2mat(M_fishrange));
        
        % check M_ClusterIDs
        ClusterIDs = M_ClusterIDs{i_set};        
        nClus = length(VAR(i_fish).ClusGroup{ClusterIDs(1)});
        if nClus>= ClusterIDs(2)
            flag2 = ~isempty(VAR(i_fish).ClusGroup{ClusterIDs(1)}.gIX);
        else
            flag2 = 0;
        end
        
        % check M_stimrange
        stimrange = M_stimrange{i_fish};
        flag3 = ~isempty(stimrange);
        
        isLoadFish = flag1&flag2&flag3;
        if ~isLoadFish
            M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish);
        end
    end
end

fishrange_load = unique(cell2mat(M_fishrange_im));
end