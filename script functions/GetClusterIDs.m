function ClusterIDs = GetClusterIDs(option)

if ~exist('option','var'),
    % default range:
    ClusterIDs = [6,1];
elseif strcmp(option,'0.7'),
    ClusterIDs = [6,1];
elseif strcmp(option,'0.5'),
    ClusterIDs = [7,1];
elseif strcmp(option,'all'),
    ClusterIDs = [2,1];
 
end

end