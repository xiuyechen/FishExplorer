% Batch run autoclusting, Sweep K2 with 2-fold CV
% Need to run Batch_k20_CV12_AK first

%% Initialize Parameters 

range_fish = 1:15; % select which fish to analyze (can be multiple)
M_stim = 1; % select which stimuli to use when clustering
k2_sweep = round(logspace(1,3,6)); % values of k2 to sweep

isFullData = 1; % If using all cells

mergeDef = 0.7; % Default values of clusParams
capDef = 0.7;
reg1Def = 0.7;
reg2Def = 0.7;
minSizeDef = 10;
k1Def = 20;
k2Def = 20;

clusParams = struct('merge',mergeDef,'cap',capDef,'reg1',reg1Def,...
    'reg2',reg2Def,'minSize',minSizeDef,'k1',k1Def,'k2',k2Def);

k2_nClus = zeros(length(k2_sweep),length(range_fish));
k2_CVscore = zeros(length(k2_sweep),length(range_fish));
k2_CVscore_raw = cell(length(k2_sweep),length(range_fish));
k2_nCells = zeros(length(k2_sweep),length(range_fish));
k2_compTime = zeros(length(k2_sweep),length(range_fish));

k2_data = struct('clusParams',clusParams,'nClus',k2_nClus,'CVscore',k2_CVscore,...
    'nCells',k2_nCells,'compTime',k2_compTime,'clusA',[],'clusB',[]);

data_masterdir = GetCurrentDataDir();

%%
disp(['Batch Script: Sweep k2, iFish']);
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
    for k2_dummy = 1:length(k2_sweep),
        
        k2start = tic;
        clusParams.k2 = k2_sweep(k2_dummy);
        Score = zeros(1,2);%(length(M_stim),2);
        %%
        nClus = zeros(1,2);
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
        k2_CVscore_raw{k2_dummy,k_fish} = Score;

        % Save clusters and scores

        k2_CVscore(k2_dummy,k_fish) = mean(Score);
        k2_nClus(k2_dummy,k_fish) = mean([nClus(1),nClus(2)]);        
        k2_nCells(k2_dummy,k_fish) = mean([nCells(1),nCells(2)]);
        k2_compTime(k2_dummy,k_fish) = toc(k2start);
        
        k2_data(k2_dummy,k_fish).CVscore = mean(Score);
        k2_data(k2_dummy,k_fish).nClus = mean([nClus(1),nClus(2)]);        
        k2_data(k2_dummy,k_fish).nCells = mean([nCells(1),nCells(2)]);
        k2_data(k2_dummy,k_fish).clusParams = clusParams;
        k2_data(k2_dummy,k_fish).compTime = toc(k2start);
        
        k2_data(k2_dummy,k_fish).clusA = struct('cIX',CIX{1},'gIX',GIX{1});
        k2_data(k2_dummy,k_fish).clusB = struct('cIX',CIX{2},'gIX',GIX{2});
        end
end
scriptTime = toc(scriptStart);
disp(['Batch Script Finished, Took ' num2str(scriptTime) ' sec']);

%%
figure;
subplot(3,1,1)
plot(k2_sweep*20,k2_nClus,'-o')
% legend('Fish8','Fish9');
%xlabel('total k for kmeans')
ylabel('# of auto-clusters')
subplot(3,1,2)
plot(k2_sweep*20,k2_CVscore,'-o')
ylim([0,1])

xlabel('total k for kmeans')
ylabel('CV (overlapping cell %)')
subplot(3,1,3)
plot(k2_sweep*20,k2_nCells,'-o');
%xlabel('total k for kameans')
ylabel('total # of cells included in clusters')

legend(leg,'location','best');