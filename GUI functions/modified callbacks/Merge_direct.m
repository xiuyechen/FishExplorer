function [gIX, numU] = Merge_direct(f,thres_merge,M_0,cIX,gIX)
M = M_0(cIX,:);
[gIX, numU] = f.HierClus(M,gIX);
U = unique(gIX);
M = M_0(cIX,:);
[C,D] = FindCentroid_Direct(gIX,M);
i = 1;
while i<numU,
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i); %#ok<*FNDSB>
        U = unique(gIX);
        numU = length(U);
        
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
        C(i+1,:) = [];
        D(i+1) = [];
    else
        i = i+1;
    end
end
[gIX, numU] = f.HierClus(M,gIX);
disp('merging complete');
end