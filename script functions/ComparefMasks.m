function [TF,prcOverlap] = ComparefMasks(fMask1,fMask2)
thres_ratio = 0.4;

IX1 = find(fMask1); IX2 = find(fMask2);
IX_int = intersect(IX1,IX2);
ratio1 = length(IX_int)/length(IX1);
ratio2 = length(IX_int)/length(IX2);

prcOverlap = max([ratio1,ratio2]);
if prcOverlap>thres_ratio,
    TF = 1;
else
    TF = 0;
end

end