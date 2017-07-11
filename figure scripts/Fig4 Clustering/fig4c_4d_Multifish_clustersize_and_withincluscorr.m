% fig2d size histogram for all clusters within one fish

% load all clusters (Fish6_woA_Master0.7)
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = 1:18;

M_Count = cell(length(range_fish),1);
M_Coeff = cell(length(range_fish),1);

for i_fish = range_fish
    ClusterIDs = [6,1];
    [cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    %% cluster sizes
    U = unique(gIX);
    numU = length(U);
    Count = zeros(numU,1);
    for i=1:numU,
        IX = find(gIX == U(i));
        Count(i) = length(IX);
    end
    M_Count{i_fish} = Count;
    
    %% within-cluster correlations
    U = unique(gIX);
    numU = length(U);
    B = zeros(numU,1);
    for i=1:numU,
        IX = find(gIX == U(i));
        coeffs = corr(M(IX,:)');
        m = coeffs(:);
        B(i) = mean(m);
    end
    M_Coeff{i_fish} = B;
    
end

%% fig4c: Distribution of cluster sizes
% make log-scaled bins for the histogram
% bins = 10:10:2100; % enough for Autoclus0.7
binlim = [10,1000];
logbins = log(binlim);
numbins = 10;
linbins = linspace(logbins(1),logbins(end),numbins+1);
xbins = exp(linbins);

% pool
N = zeros(length(range_fish),numbins);
for i_fish = range_fish
    [N(i_fish,:),edges,bin] = histcounts(M_Count{i_fish},xbins);
end

% plot
figure('Position',[500,500,300,150]);hold on; % [500,500,350,300]
for i = 1:size(N,2)
    meanN = mean(N(:,i));
    semN = std(N(:,i))/sqrt(length(range_fish));
    plot([edges(i),edges(i)],[0,meanN],'color',[1,0.5,0.5],'linewidth',10)
    plot([edges(i),edges(i)],[meanN-semN,meanN+semN],'color',[0.2,0.2,0.2],'linewidth',0.5);
end
% bar(edges(1:length(N)),N)
set(gca,'XScale','log')
xlim([7,10^3])
set(gca,'XTick',[10,100,1000])
ylim([0,50])
xlabel('cluster size')
ylabel('count')

%% fig4d: Distribution of within-cluster correlations
figure('Position',[500,500,300,150]);hold on;
xbins = 0:0.05:1;

% pool
N = zeros(length(range_fish),length(xbins)-1);
for i_fish = range_fish
    [N(i_fish,:),edges,bin] = histcounts(M_Coeff{i_fish},xbins);
end

for i = 1:size(N,2)
    meanN = mean(N(:,i));
    semN = std(N(:,i))/sqrt(length(range_fish));
    plot([edges(i),edges(i)],[0,meanN],'color',[1,0.5,0.5],'linewidth',7.5)
    plot([edges(i),edges(i)],[meanN-semN,meanN+semN],'color',[0.2,0.2,0.2],'linewidth',0.5);
end

% 
% h = findobj(gca,'Type','patch');
% h.FaceColor = [0.5 0.5 0.5];
% h.EdgeColor = 'w';
xlim([0,1])
set(gca,'XTick', 0:0.2:1);
xlabel('average corr. within cluster')
ylabel('count')
ylim([0,60])

%% fig4a... kmeans
% i_ClusGroup = 2;% 3;
% i_Cluster = 8; % 1;
% [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
% M = UpdateIndices_Manual(hfig,cIX,gIX,numU);
%
% %% left plot
% figure('Position',[50,100,800,1000]);
% % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
% setappdata(hfig,'isPlotBehavior',1);
% setappdata(hfig,'isStimAvr',0);
% UpdateTimeIndex(hfig);
% DrawTimeSeries(hfig,cIX,gIX);

%%
