% fig2d size histogram for all clusters within one fish 

% load all clusters (Fish6_woA_Master0.7)
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 6;
[cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig);

%% fig4c: Distribution of cluster sizes
U = unique(gIX);
numU = length(U);
C = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));    
    C(i) = length(IX);
end
% [N,edges] = histcounts(C,10:10:2100);
% 
% figure('Position',[500,500,350,300]);hold on;
% for i = 1:length(N),
%     plot([edges(i),edges(i)],[0,N(i)],'color',[1,0.5,0.5],'linewidth',6)
% end
% set(gca,'XScale','log')

% try making log-scaled bins for the histogram
% bins = 10:10:2100; % enough for Autoclus0.7
bins = 10:100:10700;
logbins = log(bins);
linbins = linspace(logbins(1),logbins(end),10);
xbins = exp(linbins);

[N,edges,bin] = histcounts(C,xbins);
%
figure('Position',[500,500,200,200]);hold on; % [500,500,350,300]
for i = 1:length(N),
    plot([edges(i),edges(i)],[0,N(i)],'color',[1,0.5,0.5],'linewidth',8)
end
% bar(edges(1:length(N)),N)
set(gca,'XScale','log')
xlim([5,10^4])
set(gca,'XTick',[10,100,1000,10000])
% ylim([0,70])
xlabel('cluster size')
ylabel('count')

%% fig4d: Distribution of within-cluster correlations
U = unique(gIX);
numU = length(U);
B = zeros(numU,3);
for i=1:numU,
    IX = find(gIX == U(i));
    coeffs = corr(M(IX,:)');
    m = coeffs(:);
%     B(i,1) = min(m);
    B(i,2) = mean(m);
%     B(i,3) = median(m);
end

figure('Position',[500,500,250,150]);
hist(B(:,2),0:0.05:1)
h = findobj(gca,'Type','patch');
h.FaceColor = [0.5 0.5 0.5];
h.EdgeColor = 'w';
xlim([0,1])
set(gca,'XTick', 0:0.2:1);
xlabel('average corr. within cluster')
ylabel('count')

%% fig4a... kmeans
i_ClusGroup = 2;% 3;
i_Cluster = 8; % 1;
[cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
M = UpdateIndices_Manual(hfig,cIX,gIX,numU);

%% left plot
figure('Position',[50,100,800,1000]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',1);
setappdata(hfig,'isStimAvr',0);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);

%%
