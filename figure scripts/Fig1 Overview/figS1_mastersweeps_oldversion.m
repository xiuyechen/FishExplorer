% plot masterthres sweeps
load('C:\Users\Xiu\Dropbox\FishExplorer2\AK Test Scripts\mydata.mat');

%% number of clusters
m = masterThresh_data;
data = zeros(size(m));
for i = 1:size(m,1),
    for j = 1:size(m,2),
        data(i,j) = m(i,j).nClus;
    end
end

figure('Position',[500,400,250,200]);
plot(0.5:0.05:0.9,data)
ylabel('# of clusters')
xlabel('clustering threshold (~correlation)')
xlim([0.5,0.9])
ylim([0,360])
title('custom tuning of stringency')

%% CV
data = zeros(size(m));
for i = 1:size(m,1),
    for j = 1:size(m,2),
        data(i,j) = m(i,j).CVscore;
    end
end

figure('Position',[500,400,300,200]);
plot(0.5:0.05:0.9,data)
ylabel('cross-val. score')
xlabel('clustering thresh. (~corr.)')
xlim([0.5,0.9])
ylim([0,1])
title('cross-val. (overlapping cell %)')

%% number of cells included
data = zeros(size(m));
for i = 1:size(m,1),
    for j = 1:size(m,2),
        data(i,j) = m(i,j).nCells;
    end
end

figure('Position',[500,400,300,200]);
plot(0.5:0.05:0.9,data)
ylabel('# of cells')
xlabel('clustering thresh. (~corr.)')
xlim([0.5,0.9])
title('total cell # in clusters')

