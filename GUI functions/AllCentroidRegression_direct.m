function [cIX,gIX,numK] = AllCentroidRegression_direct(M_0,thres_reg,Reg)
% M_0 = getappdata(hfig,'M_0');
% thres_reg = getappdata(hfig,'thres_reg');
% 
% Reg = FindCentroid(hfig);

Corr = corr(Reg',M_0');

[corr_max,IX] = max(Corr,[],1);
cIX = find(corr_max>thres_reg)';
gIX = IX(cIX)';
numK = length(unique(gIX));
end