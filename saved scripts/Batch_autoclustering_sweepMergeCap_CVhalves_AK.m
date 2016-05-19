% batch run full-clustering on all fish, sweep merge and cap with 2-fold CV
% Need to run Batch_k20_CV12_AK first

data_masterdir = GetCurrentDataDir();

% range_fish = [5,6,7];
% M_ClusGroup = [2,2,2,2];
% M_Cluster = [1,1,1,1];
%range_fish = 2:3;
 i_fish = 2;
% M_ClusGroup = 2;
% M_Cluster = 3;
M_stim = 1;
% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2];

%%
isFullData = 1;

mergeDef = 0.6;
capDef = 0.6;
reg1Def = 0.7;
reg2Def = 0.7;
minSizeDef = 10;
k1Def = 20;
k2Def = 50;

clusParams = struct('merge',mergeDef,'cap',capDef,'reg1',reg1Def,...
    'reg2',reg2Def,'minSize',minSizeDef,'k1',k1Def,'k2',k2Def);

mergeSweep = 0.5:.1:0.9;
capSweep = 0.5:.1:0.9;
%M_param = 2:8:10;%0.3:0.1:0.8;

Param_nClus = zeros(length(mergeSweep),length(capSweep));
Param_CVscore = zeros(length(mergeSweep),length(capSweep));
Param_CVscore_raw = cell(length(mergeSweep),length(capSweep));

%     thres_split = M_param(k_param);
%     setappdata(hfig,'thres_split',thres_split);

%for k_fish = 1:length(range_fish),
    %i_fish = range_fish(k_fish);
    leg{k_fish} = ['Fish ' num2str(range_fish(k_fish))];
    disp(i_fish);
    LoadFullFish(hfig,i_fish,isFullData);
    absIX = getappdata(hfig,'absIX');
    
    %% partitions for CV
    timelists = getappdata(hfig,'timelists');
    timelists_names = getappdata(hfig,'timelists_names');
    periods = getappdata(hfig,'periods');
 %   if length(periods)>1,
        timelistsCV = cell(length(M_stim),2);
        
        k_stim = 1;
        for k_stim = 1:length(M_stim),
            i_stim = M_stim(k_stim);
            TL = timelists{i_stim};
            period = periods(i_stim);
            nrep = size(TL,2)/periods(i_stim); % integer
            n = floor(nrep/2);
            timelistsCV{k_stim,1} = TL(1):TL(n*period);
            timelistsCV{k_stim,2} = TL(1+n*period):TL(2*n*period);
        end
 %   end
 
 for mergeDummy = 1:length(mergeSweep),
     clusParams.merge = mergeSweep(mergeDummy);
     
     for capDummy = 1:length(capSweep),
         clusParams.cap = capSweep(capDummy);
         
         Score = zeros(1,2);%(length(M_stim),2);
         %%
         NumClus = zeros(1,2);
         CIX = cell(1,2);
         GIX = cell(1,2);
         for k = 1:2, % CV halves, load them from saved
             %%
             i_ClusGroup = 2;
             i_Cluster = k+2;
             [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
             
             tIX = timelistsCV{k_stim,k};
             M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
             
             isWkmeans = 0;
             [cIX,gIX] = AutoClusteringAK(cIX,gIX,M_0,isWkmeans,clusParams);
             %%
             NumClus(k) = length(unique(gIX));
             CIX{k} = cIX;
             GIX{k} = gIX;
         end
         % plot cell-matching figure
         Score(1) = HungarianCV(NumClus(1),NumClus(2),CIX{1},CIX{2},GIX{1},GIX{2});% true,timelists_names{i_stim});
         Score(2) = HungarianCV(NumClus(2),NumClus(1),CIX{2},CIX{1},GIX{2},GIX{1});% true,timelists_names{i_stim});
         Param_CVscore(mergeDummy, capDummy) = mean(Score);
         Param_CVscore_raw{mergeDummy,capDummy} = Score;
         
         nClus1 = length(unique(GIX{1}));
         nClus2 = length(unique(GIX{2}));
         Param_nClus(mergeDummy, capDummy) = mean([nClus1,nClus2]);
     end
end
     %end

%%
figure;
[C,h] = contourf(mergeSweep,capSweep,Param_nClus);
clabel(C,h);
xlabel('Merge Threshold');
ylabel('Cap Threshold');
colorbar;
grid on;

figure;
[C2, h2] = contourf(mergeSweep,capSweep,Param_CVscore);
clabel(C2,h2);
colorbar;caxis([0 1])
xlabel('Merge Threshold');
ylabel('Cap Threshold');
grid on;
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