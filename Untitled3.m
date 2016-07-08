TF = zeros(size(Corr));
TF(Corr>thres_reg) = 1;
S = sum(TF,2);