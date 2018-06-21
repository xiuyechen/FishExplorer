% Batch run autoclusting, Sweep master threshold with 2-fold CV
% ie all clustering thresholds are the same
% Need to run Batch_k20_CV12_AK first

%% Initialize Parameters 

range_fish = 8:15; % select which fish to analyze (can be multiple)
M_stim = 2; % select which stimuli to use when clustering
masterThresh_sweep = 0.5:.1:0.9; % values of thresh to sweep

isFullData = 1; % If using all cells

masterDef = 0.7; % Default values of clusParams
mergeDef = masterDef;
capDef = masterDef;
reg1Def = masterDef;
reg2Def = masterDef;
minSizeDef = 10;
k1Def = 20;
k2Def = 500;

clusParams = struct('merge',mergeDef,'cap',capDef,'reg1',reg1Def,...
    'reg2',reg2Def,'minSize',minSizeDef,'k1',k1Def,'k2',k2Def);

masterThresh_nClus = zeros(length(masterThresh_sweep),length(range_fish));
masterThresh_CVscore = zeros(length(masterThresh_sweep),length(range_fish));
masterThresh_CVscore_raw = cell(length(masterThresh_sweep),length(range_fish));
masterThresh_nCells = zeros(length(masterThresh_sweep),length(range_fish));
masterThresh_compTime = zeros(length(masterThresh_sweep),length(range_fish));

masterThresh_data = struct('clusParams',clusParams,'nClus',masterThresh_nClus,'CVscore',masterThresh_CVscore,...
    'nCells',masterThresh_nCells,'compTime',masterThresh_compTime,'clusA',[],'clusB',[]);

data_masterdir = GetCurrentDataDir();

%%
disp(['Batch Script: Sweep masterThresh, iFish']);
scriptStart = tic;
clear leg;
for k_fish = 1:length(range_fish),
    
    %% Load fish data
    i_fish = range_fish(k_fish);
    leg{k_fish} = ['Fish ' num2str(range_fish(k_fish))];
    LoadFullFish(hfig,i_fish,isFullData); % Load the data
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
    
    %% Sweep through k2 and do clustering
    for k_masterThresh = 1:length(masterThresh_sweep),
        
        sweepStart = tic;
        clusParams.merge = masterThresh_sweep(k_masterThresh);
        clusParams.cap = masterThresh_sweep(k_masterThresh);
        clusParams.reg1 = masterThresh_sweep(k_masterThresh);
        clusParams.reg2 = masterThresh_sweep(k_masterThresh);
        Score = zeros(1,2);%(length(M_stim),2);
        %%
        nClus = zeros(1,2);
        nCells = zeros(1,2);
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
            [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams);
%%
            nClus(half) = length(unique(gIX));
            nCells(half) = length(unique(cIX));
            CIX{half} = cIX;
            GIX{half} = gIX;
        end
        
        % Perform 2x Cross Validation
        Score(1) = HungarianCV(nClus(1),nClus(2),CIX{1},CIX{2},GIX{1},GIX{2});% true,timelists_names{i_stim});
        Score(2) = HungarianCV(nClus(2),nClus(1),CIX{2},CIX{1},GIX{2},GIX{1});% true,timelists_names{i_stim});
        masterThresh_CVscore_raw{k_masterThresh,k_fish} = Score;

        % Save clusters and scores

        masterThresh_CVscore(k_masterThresh,k_fish) = mean(Score);
        masterThresh_nClus(k_masterThresh,k_fish) = mean([nClus(1),nClus(2)]);        
        masterThresh_nCells(k_masterThresh,k_fish) = mean([nCells(1),nCells(2)]);
        masterThresh_compTime(k_masterThresh,k_fish) = toc(sweepStart);
        
        masterThresh_data(k_masterThresh,k_fish).CVscore = mean(Score);
        masterThresh_data(k_masterThresh,k_fish).nClus = mean([nClus(1),nClus(2)]);        
        masterThresh_data(k_masterThresh,k_fish).nCells = mean([nCells(1),nCells(2)]);
        masterThresh_data(k_masterThresh,k_fish).clusParams = clusParams;
        masterThresh_data(k_masterThresh,k_fish).compTime = toc(sweepStart);
        
        masterThresh_data(k_masterThresh,k_fish).clusA = struct('cIX',CIX{1},'gIX',GIX{1});
        masterThresh_data(k_masterThresh,k_fish).clusB = struct('cIX',CIX{2},'gIX',GIX{2});
        end
end
scriptTime = toc(scriptStart);
disp(['Batch Script Finished, Took ' num2str(scriptTime) ' sec']);

%%
figure;
subplot(2,2,1)
plot(masterThresh_sweep,masterThresh_nClus,'-o')
% legend('Fish8','Fish9');
%xlabel('total k for kmeans')
ylabel('# of auto-clusters')

subplot(2,2,2)
plot(masterThresh_sweep,masterThresh_CVscore,'-o')
ylim([0,1])
ylabel('CV (overlapping cell %)')

subplot(2,2,3)
plot(masterThresh_sweep,masterThresh_nCells,'-o');
%xlabel('total k for kameans')
ylabel('total # of cells included in clusters')
xlabel('master threshold')

subplot(2,2,4)
plot(masterThresh_sweep,masterThresh_nCells.*masterThresh_CVscore,'-o');
%xlabel('total k for kameans')
ylabel('total # of cells passing CV')
xlabel('master threshold')
legend(leg,'location','best');