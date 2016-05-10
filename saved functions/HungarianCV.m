function score = HungarianCV(numClus1,numClus2,cIX1,cIX2,gIX1,gIX2,name)
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
figure;
imagesc(im1)
title(name);
colormap(bluewhitered)
axis equal; axis tight

score = trace(im1)/sum(sum(im1));
end