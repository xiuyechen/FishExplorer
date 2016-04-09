function [cIX,gIX,numK] = FindCellsFromCentroids(clusIX,cIX,gIX)
if size(gIX,1)<size(gIX,2),
    gIX = gIX';
end

U = unique(gIX);
range = U(clusIX)';
IX = [];
for i = range,
    IX = [IX;find(gIX==i)];
end
cIX = cIX(IX);
gIX = gIX(IX);
numK = length(clusIX);
end