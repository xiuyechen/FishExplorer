function [gIX, numK] = SqueezeGroupIX(gIX)
U = unique(gIX);
numK = length(U);
for i = 1:numK,
    old = U(i);
    gIX(gIX==old) = i;
end
end