function [C,D] = FindCentroid(hfig) % used to store C,D until reset
gIX = getappdata(hfig,'gIX');
% C = getappdata(hfig,'Centroids');
% D = getappdata(hfig,'D_ctrd');
% if isempty(C),
M = getappdata(hfig,'M');

U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
setappdata(hfig,'Centroids',C);
setappdata(hfig,'D_ctrd',D);
% end
end