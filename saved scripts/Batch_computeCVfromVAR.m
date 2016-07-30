% Batch_computeCVfromVAR.m

% range_fish = [12:15,17:18];       

%% CV within stim type        
% range_fish = [12:15,17:18]; 
% M_clus = {1,2,3,4,5,6};
%  
% nTypes = length(M_clus);
% Score = zeros(length(range_fish),nTypes,2);
% 
% for i_fishnum = 1:length(range_fish),
%     i_fish = range_fish(i_fishnum);
%     disp(i_fish);
%     
%     for i_type_ref = 1:nTypes,
%         
%         [cIX1,gIX1] = LoadCluster_Direct(i_fish,4,M_clus{i_type_ref});
%         [cIX2,gIX2] = LoadCluster_Direct(i_fish,5,M_clus{i_type_ref});
% 
%         if ~isempty(cIX1) && ~isempty(cIX2),
%             Score(i_fishnum,i_type_ref,1) = HungarianCV(cIX1,cIX2,gIX1,gIX2);%,isPlotFig,'defStim');
%             Score(i_fishnum,i_type_ref,2) = HungarianCV(cIX2,cIX1,gIX2,gIX1);
%         else
%             Score(i_fishnum,i_type_ref,1) = NaN;
%             Score(i_fishnum,i_type_ref,2) = NaN;
%         end
%     end    
% end

%%
% data_dir = GetCurrentDataDir();
% save(fullfile(data_dir,'Score_Fset3_withinStimtype.mat'),'Score');

%% CV across stim types
t = tic;
range_fish = [12:15,17:18]; 
M_clus = {1,2,3,4,5,6};

nTypes = length(M_clus);
Score_triangle = zeros(nTypes,nTypes,length(range_fish),2);

for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    disp(i_fish);
    
    for i_type_ref = 1:nTypes,
        for i_type_test = i_type_ref:nTypes, % only do half the matrix
            
            [cIX1,gIX1] = LoadCluster_Direct(i_fish,4,M_clus{i_type_ref});
            [cIX2,gIX2] = LoadCluster_Direct(i_fish,5,M_clus{i_type_test});
            
            if ~isempty(cIX1) && ~isempty(cIX2),
                score1 = HungarianCV(cIX1,cIX2,gIX1,gIX2);%,isPlotFig,'defStim');
                score2 = HungarianCV(cIX2,cIX1,gIX2,gIX1);
                Score_triangle(i_type_ref,i_type_test,i_fishnum,1) = score1;
                Score_triangle(i_type_ref,i_type_test,i_fishnum,2) = score2;
            else
                Score_triangle(i_type_ref,i_type_test,i_fishnum,1) = NaN;
                Score_triangle(i_type_ref,i_type_test,i_fishnum,2) = NaN;
            end
        end
    end
end

data_dir = GetCurrentDataDir();
save(fullfile(data_dir,'CVscores_Fset3_multistimtypes.mat'),'Score_triangle');
t = toc

%% Plot triangle matrix of CV scores between different stim types
S = mean(mean(Score_triangle,4),3);
for i=1:6,S(i,i)=Score_fishavr(i);end

M_clusnames = {'all','phototaxis','OMR','spontaneous','dark flashes','looming'};
figure;CorrPlot(S,1,M_clusnames);
set(gca,'XAxisLocation','top','XTickLabel',M_clusnames,'XTickLabelRotation',45);

% then save as svg and edit in Illustrator

%% Plot: each x is a type, each dot from a single fish
figure;
hold on;
for i_type_ref = 1:nTypes;
    
    y = mean(squeeze(Score(:,i_type_ref,:)),2); % average of the two (directional) CV's
    plot(i_type_ref*ones(size(y)),y,'o')
    ylim([0,1])
    
end

%%
Score_fishavr = mean(mean(Score,3),1);
figure;bar(Score_fishavr)

%% Manual temp
i_fish = 6;
[cIX1,gIX1] = LoadCluster_Direct(i_fish,4,1);
[cIX2,gIX2] = LoadCluster_Direct(i_fish,5,1);
[score,~,cIX_int,gIX_A,gIX_B] = HungarianCV(cIX1,cIX2,gIX1,gIX2);%,isPlotFig,'defStim');

cIX = cIX_int;
gIX = gIX_A;
%%
thres_minsize = 10;
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
