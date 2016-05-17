numNoisePoints = 40;
dots = rand(numNoisePoints,2);

numClus = 2;
clusCtr = rand(numClus,2);

stdMin = .04;
stdMax = .06;
clusStdX = rand(numClus)*stdMax+stdMin;
clusStdY = rand(numClus)*stdMax+stdMin;

sizeMin = 30;
sizeMax = 50;
clusSize = randi(sizeMax, numNoisePoints)+sizeMin;

for i = 1:numClus
    clusDotsX = normrnd(clusCtr(i,1),clusStdX(i),clusSize(i),1);
    clusDotsY = normrnd(clusCtr(i,2),clusStdY(i),clusSize(i),1);
    clusDots = [clusDotsX clusDotsY];
    dots = vertcat(dots, clusDots);
end

figure;
scatter(dots(:,1),dots(:,2),5,'MarkerFaceColor','k','MarkerEdgeColor','k');

%% kmeans 1
numK1 = 15;
idx_k1 = kmeans(dots,numK1);
colorsK1 = lines(numK1);
figure; hold on;
for i = 1:length(dots)
    scatter(dots(i,1),dots(i,2),5,'MarkerFaceColor',colorsK1(idx_k1(i),:),...
        'MarkerEdgeColor',colorsK1(idx_k1(i),:));
end

%% kmeans 2
numK2 = 5;
figure; hold on;
for j = 1:numK1
    dots_k1 = dots(idx_k1 == j,:);
    idx_k2 = kmeans(dots_k1,numK2);
    colorsK2 = lines(numK2*numK1);
    for i = 1:length(dots_k1)
        ci = (j-1)*numK1+idx_k2(i);
    scatter(dots_k1(i,1),dots_k1(i,2),5,'MarkerFaceColor',colorsK2(idx_k2(ci),:),...
        'MarkerEdgeColor',colorsK2(idx_k2(ci),:));
    end
end
    

    