function [cIX_l,gIX_l,cIX_r,gIX_r] = DivideCellsbyHemisphere(CellXYZ_norm,absIX,cIX,gIX)

% number of pixels in left/right direction for anat-yx_norm: 621
Y = CellXYZ_norm(absIX(cIX),2);
y_border = 310;
IX_l = find(Y<=y_border);
IX_r = find(Y>y_border);
cIX_l = cIX(IX_l);
cIX_r = cIX(IX_r);

if exist('gIX','var'),
    gIX_l = gIX(IX_l);
    gIX_r = gIX(IX_r);
else
    gIX_l = ones(size(cIX_l));
    gIX_r = 2*ones(size(cIX_l));
end

end