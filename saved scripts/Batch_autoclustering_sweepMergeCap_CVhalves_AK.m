% Batch run autoclusting, Sweep merge and cap with 2-fold CV
% Need to run Batch_k20_CV12_AK first

%% Initialize Parameters

i_fish = 3; % select fish number (just one)
M_stim = 1; % select which stimuli to use
mergeSweep = 0.5:.1:0.9; % values of merge (corr)
capSweep = 0.5:.1:0.9; % values of cap (corr)

isFullData = 1;

mergeDef = 0.7; % Default values of clusParams
capDef = 0.7;
reg1Def = mergeDef;
reg2Def = capDef;
minSizeDef = 10;
k1Def = 20;
k2Def = 500;

clusParams = struct('merge',mergeDef,'cap',capDef,'reg1',reg1Def,...
    'reg2',reg2Def,'minSize',minSizeDef,'k1',k1Def,'k2',k2Def);

threshSweep_nClus = zeros(length(mergeSweep),length(capSweep));
threshSweep_CVscore = zeros(length(mergeSweep),length(capSweep));
threshSweep_CVscore_raw = cell(length(mergeSweep),length(capSweep));
threshSweep_nCells = zeros(length(mergeSweep),length(capSweep));
threshSweep_compTime = zeros(length(mergeSweep),length(capSweep));

threshSweep_data = struct('clusParams',clusParams,'nClus',threshSweep_nClus,'CVscore',threshSweep_CVscore,...
    'nCells',threshSweep_nCells,'compTime',threshSweep_compTime,'clusA',[],'clusB',[]);

data_masterdir = GetCurrentDataDir();

%%
disp(['Batch Script: Sweep merge, cap, for fish ' num2str(i_fish)]);
scriptStart = tic;
clear leg;

%% Load fish data

LoadFullFish(hfig,i_fish,isFullData);
absIX = getappdata(hfig,'absIX');

%% Partitions for CV
timelists = getappdata(hfig,'timelists');
timelists_names = getappdata(hfig,'timelists_names');
periods = getappdata(hfig,'periods');
    
timelistsCV = cell(length(M_stim),2);

% Divides each stimulous in half
for k_stim = 1:length(M_stim),
    i_stim = M_stim(k_stim);
    TL = timelists{i_stim};
    period = periods(i_stim);
    nrep = size(TL,2)/periods(i_stim); % integer
    n = floor(nrep/2);
    timelistsCV{k_stim,1} = TL(1):TL(n*period);
    timelistsCV{k_stim,2} = TL(1+n*period):TL(2*n*period);
end

% Sweep Merge and Cap and do clustering
 for mergeDummy = 1:length(mergeSweep),
     clusParams.merge = mergeSweep(mergeDummy);
     clusParams.reg1 = mergeSweep(mergeDummy);
     for capDummy = 1:length(capSweep),
         clusParams.cap = capSweep(capDummy);
         clusParams.reg2 = capSweep(capDummy);
         
         Score = zeros(1,2);%(length(M_stim),2);
         %%
         NumClus = zeros(1,2);
         CIX = cell(1,2);
         GIX = cell(1,2);
         for half = 1:2, % CV halves, load them from saved
             %%
             i_ClusGroup = 2;
             i_Cluster = half+2;
             % this matches the assignments 
             % in Batch_k20_CV12_AK
             [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
             
             tIX = timelistsCV{k_stim,half};
             M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
             
             isWkmeans = 0;
             [cIX,gIX] = AutoClusteringAK(cIX,gIX,M_0,isWkmeans,clusParams);
             %%

             nClus(half) = length(unique(gIX));
             nCells(half) = length(unique(cIX));
             CIX{half} = cIX;
             GIX{half} = gIX;
         end
         
         % Perform 2x Cross Validation
         Score(1) = HungarianCV(nClus(1),nClus(2),CIX{1},CIX{2},GIX{1},GIX{2});% true,timelists_names{i_stim});
         Score(2) = HungarianCV(nClus(2),nClus(1),CIX{2},CIX{1},GIX{2},GIX{1});% true,timelists_names{i_stim});
         threshSweep_CVscore_raw{mergeDummy,capDummy} = Score;
         
         % Save clusters and scores
         threshSweep_CVscore(mergeDummy, capDummy) = mean(Score);
         threshSweep_nClus(mergeDummy, capDummy) = mean([nClus(1),nClus(2)]);
         threshSweep_nCells(mergeDummy, capDummy) = mean([nCells(1),nCells(2)]);
         threshSweep_compTime(mergeDummy, capDummy) = toc(k2start);
         
         threshSweep_data(mergeDummy, capDummy).CVscore = mean(Score);
         threshSweep_data(mergeDummy, capDummy).nClus = mean([nClus(1),nClus(2)]);
         threshSweep_data(mergeDummy, capDummy).nCells = mean([nCells(1),nCells(2)]);
         threshSweep_data(mergeDummy, capDummy).clusParams = clusParams;
         threshSweep_data(mergeDummy, capDummy).compTime = toc(k2start);
         
         threshSweep_data(mergeDummy, capDummy).clusA = struct('cIX',CIX{1},'gIX',GIX{1});
         threshSweep_data(mergeDummy, capDummy).clusB = struct('cIX',CIX{2},'gIX',GIX{2});
     end
end
     %end
scriptTime = toc(scriptStart);
disp(['Batch Script Finished, Took ' num2str(scriptTime) ' sec']);
%%
figure;

subplot(2,2,1);
[C,h] = contourf(mergeSweep,capSweep,threshSweep_nClus);
clabel(C,h);
xlabel('Merge Threshold');
ylabel('Cap Threshold');
title('Num of Clusters');
colorbar;
grid on;

subplot(2,2,2);
[C2, h2] = contourf(mergeSweep,capSweep,threshSweep_CVscore);
clabel(C2,h2);
colorbar;caxis([0 1])
xlabel('Merge Threshold');
ylabel('Cap Threshold');
title('CV score');
grid on;

subplot(2,2,3);
[C2, h2] = contourf(mergeSweep,capSweep,threshSweep_nCells);
clabel(C2,h2);
xlabel('Merge Threshold');
ylabel('Cap Threshold');
title('Num of Cells');
grid on;
colorbar;

subplot(2,2,4);
[C2, h2] = contourf(mergeSweep,capSweep,threshSweep_nCells.*threshSweep_CVscore);
clabel(C2,h2);
xlabel('Merge Threshold');
ylabel('Cap Threshold');
title('Num of Crossvalidated Cells');
grid on;
colorbar;
% figure;
% subplot(2,1,1)
% plot(M_param*20,Param_nClus)
% % legend('Fish8','Fish9');
% xlabel('total k for kmeans')
% ylabel('# of auto-clusters')
% subplot(2,1,2)
% plot(M_param*20,Param_CVscore)
% ylim([0,1])
% legend(leg);
% xlabel('total k for kmeans')
% ylabel('CV (overlapping cell %)')