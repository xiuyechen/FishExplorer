% For all clusters, correlation with nearest neighbor
C = FindCentroid(hfig);
coeffs = corr(C');

A = zeros(1,length(coeffs));
for i = 1:length(coeffs),
   coeffs(i,i) = nan;
   A(i) = max(coeffs(i,:));
end

figure; hist(A)

%% Distribution of within-cluster correlations
U = unique(gIX);
numU = length(U);
B = zeros(numU,3);
for i=1:numU,
    i
    IX = find(gIX == U(i));
    coeffs = corr(M(IX,:)');
    m = coeffs(:);
    B(i,1) = min(m);
    B(i,2) = mean(m);
    B(i,3) = median(m);
end

figure;hist(B)

%% Distribution of cluster sizes
U = unique(gIX);
numU = length(U);
C = zeros(numU,1);
for i=1:numU,
    i
    IX = find(gIX == U(i));    
    C(i) = length(IX);
end

[N,edges] = histcounts(C,10:10:2100);
%%
figure;bar(edges,[0,N]+1)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('cluster size')
ylabel('number of clusters (+1 to differentiate 0 from 1)')


%% Plot Foxel Characterization

figure;
% after second kmeans
counts = hist(gIX,1:max(gIX));
hist(counts,1:1:200);
title('Fish8: cluster sizes after 2nd kmeans');
text(20,300,'mean=10;median=4;mode=1;31 clusters>200')