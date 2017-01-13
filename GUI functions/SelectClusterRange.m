function SelectClusterRange(hfig,cIX,gIX,range)
% update indices
tempI = [];
for i = range,
    tempI = [tempI;find(gIX==i)];
end
cIX = cIX(tempI);
gIX = gIX(tempI);
UpdateIndices(hfig,cIX,gIX);
RefreshFigure(hfig);
end