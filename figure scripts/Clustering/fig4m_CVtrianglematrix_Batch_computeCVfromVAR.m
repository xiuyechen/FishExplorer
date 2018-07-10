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
tic;
range_fish = [12:15,17:18]; 
M_clus = {1,2,3,4,5,6};
switch 2
    case 1 
        M_clusgroup = [4,5]; % for CV with first/second half (chunk) of data
    case 2
        M_clusgroup = [14,15]; % for (2-fold) CV with random frames
end

nTypes = length(M_clus);
Score_triangle = nan(nTypes,nTypes,length(range_fish),2); % was zeros

for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    disp(i_fish);
    
    for i_type_ref = 1:nTypes,
        for i_type_test = i_type_ref:nTypes, % only do half the matrix
            
            [cIX1,gIX1] = LoadCluster_Direct(i_fish,M_clusgroup(1),M_clus{i_type_ref});
            [cIX2,gIX2] = LoadCluster_Direct(i_fish,M_clusgroup(2),M_clus{i_type_test});
            
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
save(fullfile(data_dir,'CVscores_multistim_randframes.mat'),'Score_triangle');
toc
%%
Score_triangle(Score_triangle==0) = nan;
%%
Score = nanmean(Score_triangle,4);
Score_std = std(Score,0,3);
Score_fishavr = nanmean(nanmean(Score,3),1);

%% Plot triangle matrix of CV scores between different stim types
S = mean(mean(Score_triangle,4),3);
for i=1:6,S(i,i)=Score_fishavr(i);end

M_clusnames = {'all','phototaxis','OMR','spontaneous','dark flashes','looming'};
figure;CorrPlot(S,1,M_clusnames);
set(gca,'XAxisLocation','top','XTickLabel',M_clusnames,'XTickLabelRotation',45);

% then save as svg and edit in Illustrator

%% Plot: each x is a type, each dot from a single fish
figure('Position',[100,100,300,200]);
hold on;
for i_fishcount = 1:length(range_fish)
    y = diag(Score(:,:,i_fishcount));
    plot(y,'o')    
end
ylim([0,1])
M_clusnames = {'all','phototaxis','OMR','spontaneous','dark flashes','looming'};
set(gca,'XTick',1:nTypes,'XTickLabel',M_clusnames,'XTickLabelRotation',45)
ylabel('HungarianCV score')

%% CV across stim types, directly count intersection

range_fish = [12:15,17:18]; 
M_clus = {1,2,3,4,5,6};
switch 1
    case 1 
        M_clusgroup = [4,5]; % for CV with first/second half (chunk) of data
    case 2
        M_clusgroup = [14,15]; % for (2-fold) CV with random frames
end

nTypes = length(M_clus);
Int_triangle = nan(nTypes,nTypes,length(range_fish),3); % was zeros

for i_fishnum = 1:length(range_fish)
    i_fish = range_fish(i_fishnum);
    disp(i_fish);
    
    for i_type_ref = 1:nTypes
        for i_type_test = i_type_ref:nTypes % only do half the matrix
            
            [cIX1,gIX1] = LoadCluster_Direct(i_fish,M_clusgroup(1),M_clus{i_type_ref});
            [cIX2,gIX2] = LoadCluster_Direct(i_fish,M_clusgroup(2),M_clus{i_type_test});
            
            if ~isempty(cIX1) && ~isempty(cIX2)
                cIX_int = intersect(cIX1,cIX2);
                
                Int_triangle(i_type_ref,i_type_test,i_fishnum,1) = length(cIX_int);
                Int_triangle(i_type_ref,i_type_test,i_fishnum,2) = length(cIX1);
                Int_triangle(i_type_ref,i_type_test,i_fishnum,3) = length(cIX2);
            end
        end
    end
end

% calculate ratio?
Score1 = Int_triangle(:,:,:,1)./Int_triangle(:,:,:,2);
Score2 = Int_triangle(:,:,:,1)./Int_triangle(:,:,:,3);
Score = (Score1+Score1)/2;%mean(cat(4,Score1,Score2),4);

Score_std = std(Score,0,3);
Score_fishavr = nanmean(nanmean(Score,3),1);

% Plot triangle matrix of CV scores between different stim types
S = mean(Score,3);
% for i=1:6,S(i,i)=Score_fishavr(i);end

M_clusnames = {'all','phototaxis','OMR','spontaneous','dark flashes','looming'};
figure;CorrPlot(S,1,M_clusnames);
set(gca,'XAxisLocation','top','XTickLabel',M_clusnames,'XTickLabelRotation',45);

%% Plot: each x is a type, each dot from a single fish
figure('Position',[100,100,300,200]);
hold on;
for i_fishcount = 1:length(range_fish)
    y = diag(Score(:,:,i_fishcount));
    plot(y,'o')    
end
ylim([0,1])
M_clusnames = {'all','phototaxis','OMR','spontaneous','dark flashes','looming'};
set(gca,'XTick',1:nTypes,'XTickLabel',M_clusnames,'XTickLabelRotation',45)
ylabel('cell overlap % for CV')