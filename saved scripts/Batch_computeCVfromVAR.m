% Batch_computeCVfromVAR.m
% CV
range_fish = GetStimRange();        
% isPlotFig = false;
        
        
nTypes = 3;
Score = zeros(length(range_fish),nTypes,2);

M_clus1 = {1,1,2};
M_clus2 = {2,3,3};

for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
    for i_type = 1:nTypes,
        
        [cIX1,gIX1] = LoadCluster_Direct(i_fish,5,M_clus1{i_type});
        [cIX2,gIX2] = LoadCluster_Direct(i_fish,5,M_clus2{i_type});

        Score(i,i_type,1) = HungarianCV(cIX1,cIX2,gIX1,gIX2);%,isPlotFig,'defStim');
        Score(i,i_type,2) = HungarianCV(cIX2,cIX1,gIX2,gIX1);

    end    
end

%% Plot: each x is a type, each dot from a single fish
figure;
hold on;
for i_type = 1:nTypes;
    
    y = mean(squeeze(Score(:,i_type,:)),2); % average of the two (directional) CV's
    plot(i_type*ones(size(y)),y,'o')
    ylim([0,1])
    
end

%%

[cIX1,gIX1] = LoadCluster_Direct(i_fish,4,5);
[cIX2,gIX2] = LoadCluster_Direct(i_fish,4,6);

Score(i,i_type,1) = HungarianCV(cIX1,cIX2,gIX1,gIX2);%,isPlotFig,'defStim');
Score(i,i_type,2) = HungarianCV(cIX2,cIX1,gIX2,gIX1);