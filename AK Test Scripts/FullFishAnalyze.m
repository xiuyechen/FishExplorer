%% Load Fish Data

% Pairwise distance for each auto-cluster
i_fish = 8;

% LoadFishDataWithoutTS(hfig,i_fish);
LoadFullFish_FULL(hfig,i_fish);

cellResp = getappdata(hfig,'CellResp');

%% Correlation Distance
tic
Dpairs = pdist(cellResp,'correlation');
toc

figure;
h = histogram(1-Dpairs);%Plot corr, not corr dist
title('Correlation Distance');


