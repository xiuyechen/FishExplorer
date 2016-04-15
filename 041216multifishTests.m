%% Q1
i_fish = 9;
clusIX = 31:40;

% C = AllCentroids{i_fish}.Centroids(clusIX,:);

% perform t-test between each pair of clusters
numClus = length(clusIX);
H = zeros(numClus,numClus);
P = zeros(numClus,numClus);
for i = 1:numClus,
    Ref_xyz = AllCentroids{i_fish}.XYZn{clusIX(i)};
    for j = 1:numClus,
        Test_xyz = AllCentroids{i_fish}.XYZn{clusIX(j)};
        [h,p] = ttest2(Ref_xyz,Test_xyz,'Vartype','unequal');
        H(i,j) = sum(h);
        P(i,j) = sum(p);
    end
end
figure;imagesc(P)

%% Q2
numBSreps1 = 1;
numBSreps2 = 1;
for i = 1:numClus,
    Ref_xyz = AllCentroids{i_fish}.XYZn{clusIX(i)};
    for j = 1:numClus,
        Test_xyz = AllCentroids{i_fish}.XYZn{clusIX(j)};
        H_reps = zeros(numBSreps1,numBSreps2);
        P_reps = zeros(numBSreps1,numBSreps2);
        for k1 = 1:numBSreps1,
            for k2 = 1:numBSreps2,
                numRand = size(Test_xyz,1);
                IX = ceil(numRand*rand(1,numRand));
                BS = Test_xyz(IX,:);
                [h,p] = ttest2(Ref_xyz,BS,'Vartype','unequal');
                H_reps(k1,k2) = sum(h);
                P_reps(k1,k2) = sum(p);
            end
        end
        H(i,j) = mean(mean(H_reps));
        P(i,j) = mean(mean(P_reps));
    end
end
figure;imagesc(P)

%% Q3 with ttest without bootstrap
i_fish1 = 8;
i_fish2 = 9;
numClus1 = size(AllCentroids{i_fish1}.Centroids,1);
numClus2 = size(AllCentroids{i_fish2}.Centroids,1);

tic
H = zeros(numClus1,numClus2);
P = zeros(numClus1,numClus2);
for i = 1:numClus1,
    Ref_xyz = AllCentroids{i_fish1}.XYZn{i};
    for j = 1:numClus2,
        Test_xyz = AllCentroids{i_fish2}.XYZn{j};
        [h,p] = ttest2(Ref_xyz,Test_xyz,'Vartype','unequal');
        H(i,j) = sum(h);
        P(i,j) = sum(p);
    end
end
toc
figure;imagesc(P)

%% Q4 from scratch
i_fish1 = 8;
i_fish2 = 9;
numClus1 = size(AllCentroids{i_fish1}.Centroids,1);
numClus2 = size(AllCentroids{i_fish2}.Centroids,1);

CostMat = zeros(numClus1,numClus2);
for i = 1:numClus1,
    Ref_xyz = AllCentroids{i_fish1}.XYZn{i};
    for j = 1:numClus2,
        Test_xyz = AllCentroids{i_fish2}.XYZn{j};
        [h,p] = ttest2(Ref_xyz,Test_xyz,'Vartype','unequal');
%         if i~=j,
            CostMat(i,j) = -sum(p); % sum(h)
%         else
%             CostMat(i,j) = 0;
%         end
    end
end

%%
Pairs = zeros(numClus1,3);
for i = 1:numClus1,
    [a,b] = max(-CostMat(i,:));
    Pairs(i,:) = [i,b,a];
end


%%
figure;
assignment1 = munkres(CostMat);
range = 1:numClus1;
IX = find(assignment1>0);
im1 = -CostMat(range(IX),assignment1(IX));
imagesc(im1)
colormap(bluewhitered)
axis equal; axis tight

%%

range_fish = [8,9];
range_clus = [{77},{32}]; % good diffuse:[{99},{45}];% sandwiched:[{70},{18}];%
figure;
DrawCellsOnAnatProj_MultipleFish(hfig,AllCentroids,range_fish,range_clus,1);



%% find best match for each cluster
% Pairs = zeros(numClus1,3);
% for i = 1:numClus1,
%     [a,b] = max(CostMat(i,:));
%     Pairs(i,:) = [i,b,a];
% end
% figure;hist(Pairs(:,3),0:0.05:1)
% IX = find(Pairs(:,3)>0.8);