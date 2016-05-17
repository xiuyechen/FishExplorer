% Plot # of auto-clusters against # of foxels

isWkmeans = 0;
range_nf = 1:5:100; % for second step, starting with k20
M_nClus = zeros(length(range_nf),1);
for i = 1:length(range_nf),% increase # of foxels,
    tic
    numK2 = range_nf(i);
    [cIX_,gIX_] = AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans,numK2);
    M_nClus(i,1) = length(unique(gIX_));
    
%     Score = TwoFoldCV(hfig,numK2); % fake
    toc;
end

figure;
plot(range_nf*20,M_nClus,'-o');

%%
range_nf = 1:5:200;
M_nClus = vertcat(M_nClus,M_nClusCopy101_200);
figure;
plot(range_nf*20,M_nClus,'-o');

