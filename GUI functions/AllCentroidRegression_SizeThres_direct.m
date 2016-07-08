function [cIX,gIX,numK] = AllCentroidRegression_SizeThres_direct(M_0,thres_reg,Reg,thres_minsize)
Corr = corr(Reg',M_0');

TF = zeros(size(Corr));
TF(Corr>thres_reg) = 1;
S = sum(TF,2);
RegIX = find(S>=thres_minsize);

S2 = sum(TF(RegIX,:),1);
cIX = find(S2>0)';
[~,IX] = max(Corr(RegIX,:),[],1);
gIX = IX(cIX)';
numK = length(unique(gIX));
end