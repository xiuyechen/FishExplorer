% fig2d size histogram for all clusters within one fish 

% load all clusters (Fish6_woA_Master0.7)

%% Distribution of within-cluster correlations
U = unique(gIX);
numU = length(U);
B = zeros(numU,3);
for i=1:numU,
    i
    IX = find(gIX == U(i));
    coeffs = corr(M(IX,:)');
    m = coeffs(:);
%     B(i,1) = min(m);
    B(i,2) = mean(m);
%     B(i,3) = median(m);
end
%%
figure('Position',[500,500,250,150]);
hist(B(:,2),0:0.02:1)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'w';
xlim([0,1])
xlabel('average corr. within cluster')
ylabel('count')

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

figure('Position',[500,500,350,300]);hold on;
for i = 1:length(N),
    plot([edges(i),edges(i)],[0,N(i)],'color',[1,0.5,0.5],'linewidth',6)
end
set(gca,'XScale','log')
% set(gca,'YScale','log')
xlim([5,2000])
ylim([0,85])
xlabel('cluster size')
ylabel('count')