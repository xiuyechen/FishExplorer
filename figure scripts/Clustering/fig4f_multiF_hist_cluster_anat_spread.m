% make histogram of average within-cluster distance for all clusters,
% pooled for all fish

clear all; close all; clc;

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range_fish = 1:18;

M_Count = cell(length(range_fish),1);

%%
% i_set = 1;
for i_fish = range_fish
    %% load fish
    ClusterIDs = [6,1];
    [cIX_load,gIX_load] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    cIX_abs = absIX(cIX_load);
    M_xyz_norm = CellXYZ_norm(cIX_abs,:);               

    %% compute average within-cluster-distance for each cluster    
    numClus = length(unique(gIX_load));
    D = zeros(numClus,1);
    for i_clus = 1:numClus
        IX = find(gIX_load == i_clus);
        XYZ_clus = M_xyz_norm(IX,:);
       D(i_clus) = mean(pdist(XYZ_clus,'euclidean'));
    end

    M_Count{i_fish} = D;
end

%% fig4f: Distribution of within-cluster distance
figure('Position',[500,500,150,120]);hold on;
xbins = 0:50:400;

% pool
N = zeros(length(range_fish),length(xbins)-1);
for i_fish = range_fish
    [N(i_fish,:),edges,bin] = histcounts(M_Count{i_fish},xbins);
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
xlim([-25,max(xbins)])
% set(gca,'XTick', 0:0.2:1);
xlabel('avr within-clus dist')
ylabel('count')
ylim([0,60])