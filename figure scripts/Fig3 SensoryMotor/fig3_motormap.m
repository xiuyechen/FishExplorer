% find top 100 cells correlating to motor, screen with Rh4+Rh5 masks

% corr sweep (0.5-0.7) to get ~2000 cells, then kmeans, save; rank by stim-lock, save
%% folder setup
outputDir = GetOutputDataDir;
saveDir1 = fullfile(outputDir,'motor_map');
saveDir2 = fullfile(outputDir,'motor_map_stimlock');
if ~exist(saveDir1, 'dir'), mkdir(saveDir1), end;
if ~exist(saveDir2, 'dir'), mkdir(saveDir2), end;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
% i_fish = 8;

%% run fish
range_fish = GetFishRange;%[1:3,5:18];
M_thres_reg = zeros(1,18);
M_numTopCorr = zeros(1,18);

for i_fish = range_fish
    ClusterIDs = [1,1];
    [~,~,~,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    numTopCorr = 100;
    regressors = GetMotorRegressor(behavior,i_fish);
    
    %%
    M_lr = [1,3];
    cIX_seed = [];
    gIX_seed = [];
    for i_lr = 1:2,
        % left:
        regressor = regressors(M_lr(i_lr)).im;
        
        %% for each cell, find correlation coeff
        M_corr = corr(regressor',M_0');%cell_resp_ave');
        
        [~,I] = sort(M_corr,'descend');
        cIX = I(1:numTopCorr)';
        gIX = i_lr*ones(size(cIX));
        
        %% screen with anat masks #222 and #233 (Rhombomore 4 and 5)
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        
        Msk_IDs = [222,223]; % manual input
        [cIX,gIX] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
        
        num = numTopCorr;
        while length(cIX)<3, % repeat with larger numTopCorr
            num = num+100;
            cIX = I(1:num)';
            gIX = i_lr*ones(size(cIX));
            [cIX,gIX] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
            M_numTopCorr(i_fish) = num;
        end
        
        cIX_seed = [cIX_seed;cIX];
        gIX_seed = [gIX_seed;gIX];
    end
    M = UpdateIndices_Manual(hfig,cIX_seed,gIX_seed);
    Reg = FindClustermeans(gIX_seed,M);
    
%     figure('Position',[50,100,800,1000]);
%     I = LoadCurrentFishForAnatPlot(hfig,cIX_seed,gIX_seed);
%     DrawCellsOnAnat(I);

    %% regression, dynamic thresholding
    M_target_numcell= [1800,2000,2200];
    for i_numcell = 2,%1:3,
        
        target_numcell = M_target_numcell(i_numcell);
        
        Corr = corr(Reg',M_0');
        [corr_max,IX] = max(Corr,[],1);
        [~,I] = sort(corr_max,'descend');
        cIX = I(1:target_numcell)';
        gIX = IX(cIX)';
        M = UpdateIndices_Manual(hfig,cIX,gIX);
        
        M_thres_reg(i_fish) = corr_max(I(target_numcell));
        
        %% k-means
        numK = 16; % manual
        gIX = Kmeans_Direct(M,numK);
        
        %% Plot anat
        figure('Position',[50,100,800,1000]);
        I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
        DrawCellsOnAnat(I);
        
        %% Save plot
        filename = fullfile(saveDir1, ['Fish',num2str(i_fish),'_motormap_',num2str(target_numcell)]);
        saveas(gcf, filename, 'png');
        close(gcf)
        
        %%
        %% Rank by stim-lock
        [gIX,rankscore] = RankByStimLock_Direct(hfig,gIX);
        
        setappdata(hfig,'rankscore',round(rankscore*100)/100);
        setappdata(hfig,'clrmap_name','jet');
        
        % Plot anat
        figure('Position',[50,100,800,1000]);
        I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
        DrawCellsOnAnat(I);
        
        % Save plot
        filename = fullfile(saveDir2, ['Fish',num2str(i_fish),'_motormap_stimlock_anat']);
        saveas(gcf, filename, 'png');
        close(gcf)
        
        % Plot functional traces
        figure('Position',[50,100,800,1000]);
        % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
        setappdata(hfig,'isPlotBehavior',1);
        setappdata(hfig,'isStimAvr',0);
        UpdateTimeIndex(hfig);
        DrawTimeSeries(hfig,cIX,gIX);
        
        % Save plot
        filename = fullfile(saveDir2, ['Fish',num2str(i_fish),'_motormap_stimlock_func']);
        saveas(gcf, filename, 'png');
        close(gcf)
    end  
end

