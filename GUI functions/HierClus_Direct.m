function gIX = HierClus_Direct(C,gIX)
[gIX, numU] = SqueezeGroupIX(gIX);
% [C,~] = FindCentroid_Direct(gIX,M);
D = pdist(C,'correlation');
if size(C,1)>1,
    tree = linkage(C,'average','correlation');
    leafOrder = optimalleaforder(tree,D);
    
    % sort for uniform colorscale
    temp = zeros(size(gIX));
    for i = 1:numU,
        temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
    end
    gIX = temp;
end
end