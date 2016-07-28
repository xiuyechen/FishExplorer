% fig2i HungarianCV matrix plot

i_fish = 6;
[cIX1,gIX1] = LoadCluster_Direct(i_fish,4,5);
[cIX2,gIX2] = LoadCluster_Direct(i_fish,4,6);

% [score,im1] = HungarianCV(cIX1,cIX2,gIX1,gIX2,isPlotFig);
[score,im1] = HungarianCV(cIX2,cIX1,gIX2,gIX1);

%%
figure;
% subplot(1,2,1)
% imagesc(-im1)
% colormap(bluewhitered)
% axis equal; axis tight;axis xy
% 
% subplot(1,2,2)
imagesc(-log(im1))
% colormap(bluewhitered)
axis equal; axis tight;%axis xy
colormap('gray')
ylabel('clusters (1st half)')
xlabel('clusters (2nd half)')
title('Cross-Val: # of cells overlap')
text(size(im1,1)/4,size(im1,2)/20,'score = 0.69');

