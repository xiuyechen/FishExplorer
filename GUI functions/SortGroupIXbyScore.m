function [gIX,B] = SortGroupIXbyScore(H,gIX,numU,descend) %#ok<INUSD> 
% for ranking functions % new gIX is sorted based on H, size(H)=[numU,1];
if exist('descend','var'),
    [B,I] = sort(H,'descend');
else
    [B,I] = sort(H);
end
gIX_last = gIX;
for i = 1:numU,
    gIX(gIX_last==I(i)) = i;
end
end