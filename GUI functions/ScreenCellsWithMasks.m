function [cIX,gIX] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX)
% convert cell locations to pixel ID
X = CellXYZ_norm(absIX(cIX),1);
Y = CellXYZ_norm(absIX(cIX),2);
Z = CellXYZ_norm(absIX(cIX),3);
pxID = sub2ind([MASKs.height,MASKs.width,MASKs.Zs],X',Y',Z');

% screen cells
px_hist = full(sum(MASKs.MaskDatabase(:,Msk_IDs),2));
IX = ismember(pxID,find(px_hist));
cIX = cIX(IX);
gIX = gIX(IX);
end