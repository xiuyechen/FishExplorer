function [cIX,gIX,numK] = AllCentroidRegression_SizeThres_direct(M,cIX_reg,thres_reg,Reg,thres_minsize)
Corr = corr(Reg',M');

TF = zeros(size(Corr));
TF(Corr>thres_reg) = 1;
S = sum(TF,2);
RegIX = find(S>=thres_minsize);

S2 = sum(TF(RegIX,:),1);
cIX_lc = find(S2>0)'; % local index
[~,IX] = max(Corr(RegIX,:),[],1);
gIX = IX(cIX_lc)';
numK = length(unique(gIX));

cIX = cIX_reg(cIX_lc);
end