% fig2f
clear all;
% manual:
% load all clus (Fish6_woA_Master0.7);
% union with all valid cells;
% export to workspace
% save as .mat
load('tSNE_input_F6M0.5.mat');

%%
addpath(genpath('C:\Users\Xiu\Dropbox\t-SNE'));

tic
mappedA = tsne(M, [], 2); % no_dims = 2
toc

save('mappedA_tsne_AH.mat','mappedA');

%%
figure;
hold on;
temp_gIX = ones(size(gIX));
cmap = [0.5,0.5,0.5];
gscatter(mappedA(:,1),mappedA(:,2),gIX,cmap,'.',[],'off');

IX = find(gIX ~= max(gIX));
n = round((numK-1)*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(IX,1),mappedA(IX,2),gIX(IX),cmap,'.',[],'off');

axis off

%%
figure;
% hold on;
% temp_gIX = ones(size(gIX));
% cmap = [0.5,0.5,0.5];
% gscatter(mappedA(:,1),mappedA(:,2),gIX,cmap,'.',[],'off');
% 
% IX = find(gIX ~= max(gIX));
% n = round((numK-1)*1.1);
% cmap = hsv(max(1,n));
cmap = flipud(jet(numK));
gscatter(mappedA(IX,1),mappedA(IX,2),gIX(IX),cmap,'.',[],'off');

axis off