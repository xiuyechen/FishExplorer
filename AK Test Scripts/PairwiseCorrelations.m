% AK 20160415 Clustering Play

%clearvars -except M_0_fluo;

% Export M0_fluo all cells, around 40k
f_cells = M_0_fluo;
num_cells = length(f_cells);

num_select = 2000;
clear idx_select;
idx_select = randperm(num_cells, num_select);
f_cells_select = f_cells(idx_select,:);


%% Correlation Distance
tic
Dpairs = pdist(f_cells_select,'correlation');
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
idx = find(h.BinEdges < 0.1 & h.BinEdges > -0.1);
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

%% Cumulative Significant Correlations

fracCorrSig = (allCorrVals'-f(allBins));
%figure;plot(allBins,fracCorrSig);
cumCorrSig = cumsum(fracCorrSig);
figure;plot(allBins,cumCorrSig/cumCorrSig(end));
grid on;
