function [C,ClusRes] = FindClustermeans(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
for i = 1:numU
    IX = find(gIX == U(i));
    if length(IX)==1
        C(i,:) = M(IX,:);
    else
        M_s = M(IX,:);
        C(i,:) = mean(M_s);
    end
end

if nargout>1
    ClusRes = zeros(size(M));
    for i = 1:numU
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        ClusRes(IX,:) = M_s-repmat(C(i,:),length(IX),1);
    end
end
end