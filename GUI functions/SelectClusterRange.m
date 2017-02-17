function [cIX,gIX,IX] = SelectClusterRange(cIX,gIX,range)
IX = ismember(gIX,range);
cIX = cIX(IX);
gIX = gIX(IX);
end