function C = FindClustermeans(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
for i = 1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        C1 = mean(M_s);
        C(i,:) = C1;
    end
end
end