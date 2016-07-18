function score = HungarianCV(cIX1,cIX2,gIX1,gIX2,isPlotFig,name)
numClus1 = length(unique(gIX1));
numClus2 = length(unique(gIX2));

% CV by matching cell ID's
CostMat = zeros(numClus1,numClus2);
for i = 1:numClus1,
    for j = 1:numClus2,
        A = cIX1(gIX1==i);
        B = cIX2(gIX2==j);
        CostMat(i,j) = -length(intersect(A,B));
    end
end

assignment1 = munkres(CostMat);
range = 1:numClus1;
IX = find(assignment1>0);
im1 = -CostMat(range(IX),assignment1(IX));

if exist('isPlotFig','var'),
    if isPlotFig,
        figure;
        imagesc(im1)
        title(name);
        colormap(bluewhitered)
        axis equal; axis tight
    end
end

score = trace(im1)/sum(sum(im1));
end