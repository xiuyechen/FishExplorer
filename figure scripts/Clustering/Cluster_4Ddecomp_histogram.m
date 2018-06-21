% This script loads the AutoClus from each fish, and compute the
% 4D-components for each cluster. Summarize results in some histogram.

% Question remains whetherhow summary should be weighed by cluster size.

%%
clear all; close all; clc

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = 1:18;%GetFishRange;%[1:3,5:18];%
% ClusterIDs = [2,1];
ClusterIDs = [6,1];

Betas = cell(2,18);


%%
figure('Position',[50,50,500,300]);

for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    %% load data for chosen stim range
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault...
        (i_fish,hfig,ClusterIDs);

    %% cluster based
    gIX = gIX_load;
    C = FindClustermeans(gIX,M);
    Data = C;

    tic
    for i_lr = 1:2        
        [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
        %     [Data_tAvr,Data_tRes,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
        [motor_tAvr,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace(hfig,behavior);
        
        b1 = corr(motor_tRes(i_lr,:)',Data_p');
        b2 = corr(motor_tAvr(i_lr,:)',Data_p');
        b3 = sqrt(abs(var(Data_tAvr')./var(Data_p') - b2.^2)); % var(Data') is all 1
        
        % % sum(b1^2+b2^2+b3^2+b4^2) = 1
        b4 = sqrt(1-b1.^2-b2.^2-b3.^2);
        
        % assert: check orthogonality
        %     dot(motor_tAvr,motor_tRes,2)
        %     figure;histogram(dot(Data_tAvr,Data_tRes,2))
        
        % save
        Betas{i_lr,i_fish} = vertcat(b1,b2,b3,b4)';
        
    end
    toc
    
    %% pool for multi-fish histogram
    [n_pass{i_lr},edges] = histcounts(y_pass_n,xbins);
    
    %% cumulative plot?
    [B,IX_sort] = sort(b4,'descend');
    % get number of cells for each cluster
    Ncell_clus = histcounts(gIX_load,length(unique(gIX_load)));
    N_sorted = Ncell_clus(IX_sort);
    Y = cumsum(N_sorted);
    X = 1-B.^2;
    
    subplot(131);hold on; plot(X,Y)
    subplot(132);hold on; plot(X)
    subplot(133);hold on; plot(Y)
end
%%
subplot(131);hold on; xlabel('b1^2+b2^2+b3^2'); ylabel('cum. # of cells')
subplot(132);hold on; xlabel('# of clusters'); ylabel('b1^2+b2^2+b3^2')
subplot(133);hold on; xlabel('# of clusters'); ylabel('cum. # of cells')

%% histogram
nbins = 20;
figure; hold on
subplot(411);
h1 = histogram(b1,nbins);
subplot(412);
h2 = histogram(b2,nbins);
subplot(413);
h3 = histogram(b3,nbins);
subplot(414);
h4 = histogram(b4,nbins);

h1.FaceColor = 'g';
h2.FaceColor = 'r';
h3.FaceColor = 'b';
h4.FaceColor = 'y';

% h1 = histogram('BinEdges',edges,'BinCounts',mean(N_pass,1));
%     h2 = histogram('BinEdges',edges,'BinCounts',mean(N_int,1));
% 
%     h1.FaceColor = clr_XY(i_ax,:);
%     h1.EdgeColor = 'w';
%     h1.FaceAlpha = 0.6;
%     h2.FaceColor = [1,0.5,0.5];
%     h2.EdgeColor = 'w';
%     h2.FaceAlpha = 0.6;
%     
%     xlim([0,1])
%     xlabel('% beta')
%     ylabel('a.u.')
%     set(gca,'YTick',[]);

%% with GUI: manual select
[B,IX_sort] = sort(b4,'descend');
range_gIX = IX_sort(1:50);
[cIX,gIX] = SelectClusterRange(cIX_load,gIX_load,range_gIX);

% %% cumulative plot?
% % get number of cells for each cluster
% Ncell_clus = histcounts(gIX_load,length(unique(gIX_load)));
% N_sorted = Ncell_clus(IX_sort);
% Y = cumsum(N_sorted);
% X = 1-B.^2;
% 
% subplot(121);plot(X,Y)
% subplot(122);plot(Y)
figure;plot(X)
%% save 4D-decomposition
% save(fullfile(outputDir,'4D_SM_A0.7_betas.mat'),'Betas');
