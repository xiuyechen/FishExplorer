% AK 20160415 Clustering Play

%clearvars -except M_0_fluo;

% Export M (M0_fluo?) all cells, around 40k
%f_cells = M_0_fluo;
f_cells = M;
num_cells = length(f_cells);

num_select = 1000;%num_cells;
clear idx_select;
idx_select = randperm(num_cells, num_select);
f_cells_select = f_cells(idx_select,:);


%% Correlation Distance
tic
Dpairs = pdist(f_cells_select,'correlation');
% 138 sec for 86k cells
% 1.8 sec for 10k cells
toc
figure;
h = histogram(1-Dpairs);%Plot corr, not corr dist
title('Correlation Distance');

%% Euclidean Distance

% tic
% Dpairs = pdist(f_cells_select,'euclidean');
% toc
% 
% figure;
% histogram(Dpairs);
% title('Euclidean Distance');

%% Fit random variance
% Fit Gaussian to part above 1
%h = histogram(1-Dpairs);
idx = find(h.BinEdges < 0.03 & h.BinEdges > -0.03);
allBins = h.BinEdges(1:end-1) + h.BinWidth/2;
allCorrVals = h.Values;
unCorrBins = h.BinEdges(idx) + h.BinWidth/2;
unCorrVals = h.Values(idx+1);

figure;
options = fitoptions('gauss1', 'Lower', [-Inf 1 -Inf], 'Upper', [Inf 1 Inf]);
% Forcing b1 = 1 doesn't work well
%f = fit(unCorrBins',unCorrVals','gauss1',options);
f = fit(unCorrBins',unCorrVals','gauss1');
plot(f,allBins,allCorrVals);
%plot(f,unCorrBins,unCorrVals);
hold on;
grid on;
title('Pairwise Cell Correlations');
xlabel('Correlation')
ylabel('Frequency')

xlim([-0.2,0.2]);

%% Percentile Thresh
T = size(f_cells_select,2);
clear threshCorr T_trunc
T_trunc = 320:320:2880;
for i = 1:length(T_trunc)
    f_trunc = f_cells(:,1:T_trunc(i));
    
    num_select = 1000;%num_cells;
    clear idx_select;
    idx_select = randperm(num_cells, num_select);


% Correlation Distance
tic
Dpairs = pdist(f_trunc,'correlation');
% 138 sec for 86k cells
% 1.8 sec for 10k cells
toc
%figure;
%h = histogram(1-Dpairs);%Plot corr, not corr dist
title('Correlation Distance');
threshPerc = .99;
sortCorr = sort(1-Dpairs);
threshIDX = round(length(sortCorr)*threshPerc);
threshCorr(i) = sortCorr(threshIDX);
end

figure;plot(T_trunc,threshCorr,'o-')

%% Cumulative Correlations

fracCorr = (allCorrVals');
%figure;plot(allBins,fracCorrSig);
cumCorr = cumsum(fracCorr);
figure;plot(allBins,cumCorr/cumCorr(end));
grid on;
ylabel('Cumulative Frac of Pairwise Corr Distances');
xlabel('Correlation');
xlim([-0.2,0.2]);
%% Cumulative Significant Correlations

fracCorrSig = (allCorrVals'-f(allBins));
%figure;plot(allBins,fracCorrSig);
cumCorrSig = cumsum(fracCorrSig);
figure;plot(allBins,cumCorrSig/cumCorrSig(end));
grid on;
 %% How many isolated cells?
% 
% % tic
% % Dpairs = pdist2(f_cells_select,f_cells_select,'correlation');
% % toc
% 
% %maxNumLoad = 10000;
% tic
% DpairsSquare = squareform(Dpairs);
% sortD = sort(DpairsSquare);
% corrThresh = 1.3;
% threshD = sortD;
% threshD(sortD < corrThresh) = 0;
% threshD(sortD >= corrThresh) = 1; %there is probably a better way to do this
% toc % 3 sec for 10k cells
% 
% figure;imagesc(threshD);
% degree = sum(threshD,1);
% figure; hist(degree);
% 
% min_neighbors = 1;
% frac_iso = sum(degree < min_neighbors)/length(degree)

%% How many isolated cells (lower memory usage)

%Try again
corrThresh = 0.2;
distThresh = 1-corrThresh;
num_cells = length(f_cells_select);

tic

degree = zeros(num_cells,1);
m = num_cells;

for i = 1:num_cells
    j_less = 1:i-1;
    j_more = i+1:num_cells;
    %dists = pdist2(f_cells_select(i,:),f_cells_select,'correlation');
    D_idxs_less = (j_less-1).*(m-j_less/2)+i-j_less;
    D_idxs_more = (i-1).*(m-i/2)+j_more-i;
    D_idxs = [D_idxs_less D_idxs_more];
    dists = Dpairs(D_idxs);
    % D ij is D((i–1)*(m–i/2)+j–i) for i < j
    adjacency = zeros(1,num_cells);
    adjacency(dists <= distThresh) = 1;
    degree(i)=sum(adjacency);
end
toc%240 sec for 80k cells / 2 sec for 10k cells

%% Plot frac_iso as function of min neighbors

for min_neighbors =1:20
frac_iso(min_neighbors) = sum(degree < min_neighbors)/length(degree);
end
figure; bar(frac_iso);
xlabel('min neighbors');
ylabel('perc isolated cells');
title(strcat('corrThresh = ',num2str(corrThresh)));
axis([0 21 0 .8]);
%% Hierarchical Clustering
tic
Y = pdist(f_cells_select,'correlation');
Z = linkage(Y);
leafOrder = optimalleaforder(Z,Y);
%figure; dendrogram(Z,0,'ColorThreshold',1.2,'Reorder',leafOrder)
T = cluster(Z,'cutoff',2);
%figure;hist(Z(:,3));
toc

%% Density Clustering    
tic
corrThresh = 0.3:0.01:1;
for i = 1:length(corrThresh)
    epsilon=1-corrThresh(i);
    MinPts=3;
    
    clear IDX;

    IDX=DBSCAN_AKmod(f_cells_select,epsilon,MinPts);

    numClus(i) = length(unique(IDX))-1;
%     clear clusSize
%     for clusNum = 0:numClus-1
%         clusSize(clusNum+1) = length(find(IDX == clusNum));
%     end
%     clusSize
end
    toc
    figure;
plot(corrThresh,numClus);