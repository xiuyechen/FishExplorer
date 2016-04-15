% Cross-Validation



assignment1 = munkres(CostMat);
range = 1:numClus1;
IX = find(assignment1>0);
im1 = -CostMat(range(IX),assignment1(IX));
imagesc(im1)
colormap(bluewhitered)
axis equal; axis tight
