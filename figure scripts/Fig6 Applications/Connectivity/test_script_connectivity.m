% BCT
clear all; close all; clc

% BCT_dir = 'C:\Users\Xiu\Dropbox (Personal)\BCT';
% addpath(genpath(BCT_dir));
data_masterdir = GetCurrentDataDir();

%% load data from one fish
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

i_fish = 6;
ClusterIDs = [6,1];
[cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);

absIX = getappdata(hfig,'absIX');

%%

data = M;
%%
clusMeans = FindCentroid_Direct(gIX,M);
data = clusMeans;

%% Convert to a weighted adjacency matrix
% correlation matrix of cluster means
tic
W_full = corrcov(data*data');
toc

figure;imagesc(W_full);

% threshold and remove diagonal
% thres = 0.2;
% % % W = threshold_absolute(corrM,thres);
% % W(1:size(W,1)+1:end)=0;                 %clear diagonal
% W(W<thres)=0; 

% normalize to range 0 to 1
% W = W./max(abs(W(:)));
% W = weight_conversion(W, 'autofix');

%%
tic
D = pdist(data,'correlation');
tree = linkage(data,'average','correlation');
leafOrder = optimalleaforder(tree,D);
toc
%     
data2 = data(leafOrder,:);
%%
W_full_sorted = corrcov(data2*data2');
figure;imagesc(W_full_sorted);
axis image

% save(fullfile(data_masterdir,'graph_fish6_Auto07.mat'),'M','W_full','leafOrder','W_full_sorted');


%% threshold again
thres = 0.3;
W = W_full_sorted;
W(W<thres)=0; 
figure;imagesc(W);colormap(jet);colorbar;

%% Get anatomical locations of cluster means

[CellXYZ_norm,IX_inval_norm] = GetNormCellCord(i_fish);
selectXYZ_norm = CellXYZ_norm(absIX,:);
clusXYZ = NaN(max(gIX),3);
for i = 1:max(gIX)
    clusXYZ(i,:) = mean(selectXYZ_norm(cIX(gIX == i),:),1);
end

%% visualize 
% use anatomical location of cluster means 
% generate list for plotting edges
% [X,Y,Z] = adjacency_plot_und(A_nrm,clusXYZ);
% figure; plot3(X,Y,Z,'o-');

G = graph(W~=0,'OmitSelfLoops');
%%
figure;
p = plot(G);
%%
% color nodes by degree
G.Nodes.NodeColors = degree(G);
p.NodeCData = G.Nodes.NodeColors;
colorbar

%% or color node matching GUI
nClus = length(unique(gIX));
cmap = jet(nClus);
p.NodeColor = flipud(cmap);

%% node size by cluster size
nClus = length(unique(gIX));
clusSize = histcounts(gIX,0.5:1:nClus+0.5);
p.MarkerSize = sqrt(clusSize);

%% edge thickness by weight
maxLWidth = 4; % 7 in Matlab official example
E = G.Edges.EndNodes;
nNodes = size(G.Nodes,1);
IX = (sub2ind([nNodes,nNodes],E(:,1)',E(:,2)'))';
G.Edges.Weight = W(IX);
weights = W(IX);
G.Edges.LWidths = maxLWidth*(weights-min(weights))/(max(weights)-min(weights))+0.1;
p.LineWidth = G.Edges.LWidths;

%%
thres_degree = 10;%2000;
IX = find(degree(G)>thres_degree);
figure
H = subgraph(G,IX);
p1 = plot(H);
p1.NodeCData = H.Nodes.NodeColors;
set(p1,'MarkerSize',10);
set(p1,'EdgeColor',[0.5,0.5,0.5]);
colorbar
% p1 = plot(H,'NodeCData',H.Nodes.NodeColors,'LineWidth',H.Edges.LWidths);
%%

% edge thickness by weight
maxLWidth = 5; % 7 in Matlab official example
E = H.Edges.EndNodes;
nNodes = size(H.Nodes,1);
IX = (sub2ind([nNodes,nNodes],E(:,1)',E(:,2)'))';
H.Edges.Weight = W(IX);
H.Edges.LWidths = maxLWidth*H.Edges.Weight/max(H.Edges.Weight)+1;
p1.LineWidth = H.Edges.LWidths;

% display
set(p1,'MarkerSize',10);
set(p1,'EdgeColor',[0.5,0.5,0.5]);
labels = cell(1,nNodes);
for i=1:nNodes,labels{i} = num2str(i);end
set(p1,'NodeLabel',labels);
%%

% edge thickness by weight
% maxLWidth = 5; % 7 in Matlab official example
% E = G.Edges.EndNodes;
% nNodes = size(G.Nodes,1);
% IX = (sub2ind([nNodes,nNodes],E(:,1)',E(:,2)'))';
% G.Edges.Weight = W(IX);
% G.Edges.LWidths = maxLWidth*G.Edges.Weight/max(G.Edges.Weight);
% p.LineWidth = G.Edges.LWidths;
% 
% % display
% set(p,'MarkerSize',10);
% set(p,'EdgeColor',[0.5,0.5,0.5]);
% labels = cell(1,nNodes);
% for i=1:nNodes,labels{i} = num2str(i);end
% set(p,'NodeLabel',labels);

%%
figure
H = subgraph(G,[1:31 36:41]);
p1 = plot(H,'NodeCData',H.Nodes.NodeColors,'LineWidth',H.Edges.LWidths);
colorbar
path = shortestpath(H,13,9);
highlight(p1,[11 37])
highlight(p1,path,'EdgeColor','r')
%%

% Load reference brain outlines
load(fullfile(data_masterdir,'FishOutline.mat'));
figure;hold on;
imagesc((~outline_XY)); colormap(gray); axis equal;
plot(Y,X,'o-');

%
load(fullfile(data_masterdir,'FishOutline.mat'));
figure;hold on;
imagesc((~outline_YZ)); colormap(gray); axis equal;
plot(Z,X,'o-');

%
load(fullfile(data_masterdir,'FishOutline.mat'));
figure;hold on;
imagesc((~outline_ZX)); colormap(gray); axis equal;
plot(Y,Z,'o-');