function [gIX,cIX,IX] = SetGroupIXorder(gIX,order,cIX)
% 'order' can be a subset of gIX values

gIX0 = gIX;
gIX = zeros(size(gIX0));
for i = 1:length(order)
    ix = find(gIX0==order(i));
    gIX(ix) = i;
end
% delete clusters not specified in 'order'
IX = find(gIX~=0);
gIX = gIX(IX);
cIX = cIX(IX);

end