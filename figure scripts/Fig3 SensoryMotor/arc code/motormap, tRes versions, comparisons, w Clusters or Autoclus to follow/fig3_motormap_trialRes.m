% new procedure: generate motor seed/map in trialRes space
% find top 100 cells correlating to motor, screen with Rh4+Rh5(+Rh6??) masks

% corr sweep (0.5-0.7) to get top %2 cells, then kmeans, save; rank by stim-lock, save
clear all; close all; clc

%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
saveDir1 = fullfile(outputDir,'motor_map_030417');
saveDir2 = fullfile(outputDir,'motor_map_seed_030417');
if ~exist(saveDir1, 'dir'), mkdir(saveDir1), end;
if ~exist(saveDir2, 'dir'), mkdir(saveDir2), end;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',0);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];
M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

for i_fish = 1:18 %range_fish
    ClusterIDs = [1,1];
    [~,~,~,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    numTopCorr = 100;
    [~,~,regressor_m_raw] = GetMotorRegressor(behavior,i_fish);
    
    %% convert to trialRes (converting locally, not updating hfig to isTrialRes)
    [~,regressor_m_raw_tRes] = GetTrialAvrLongTrace(hfig,regressor_m_raw);
    [~,M_0_tRes] = GetTrialAvrLongTrace(hfig,M_0);
    
    %%
    M_lr = [1,3];
%     M_name = {'Left','Right'};
    cIX_seed = [];
    gIX_seed = [];
    for i_lr = 1:2,
        % left:
        reg_tRes = regressor_m_raw_tRes(M_lr(i_lr),:);
        
        %% for each cell, find correlation coeff
        M_corr = corr(reg_tRes',M_0_tRes');%cell_resp_ave');
        
        [~,I] = sort(M_corr,'descend');
        cIX = I(1:numTopCorr)';
        gIX = i_lr*ones(size(cIX));
        
        %% screen with anat masks #222 and #233 (Rhombomore 4 and 5)
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        
        Msk_IDs = [222,223,224]; % manual input
        [cIX,gIX] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
        
        % enforce that seed cells for the left are in the left hemisphere,
        % and vice versa for the right
        if i_lr ==1,
            [cIX,gIX]  = DivideCellsbyHemisphere(CellXYZ_norm,absIX,cIX,gIX);
        else
            [~,~,cIX,gIX] = DivideCellsbyHemisphere(CellXYZ_norm,absIX,cIX,gIX);
        end
        
        num = numTopCorr;
        % need at least 10 cells to qualify as a seed
        while length(cIX)<10, % repeat with larger numTopCorr
            num = num+100;
            cIX = I(1:num)';
            gIX = i_lr*ones(size(cIX));
            [cIX,gIX] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
            M_numTopCorr(i_fish) = num;
            
            % screen by left/right again
            if i_lr ==1,
                [cIX,gIX]  = DivideCellsbyHemisphere(CellXYZ_norm,absIX,cIX,gIX);
            else
                [~,~,cIX,gIX] = DivideCellsbyHemisphere(CellXYZ_norm,absIX,cIX,gIX);
            end
        end
        
        cIX_seed = [cIX_seed;cIX];
        gIX_seed = [gIX_seed;gIX];
        
    end
    % this is not tRes, just loading default:
    M = UpdateIndices_Manual(hfig,cIX_seed,gIX_seed);
    Reg = FindClustermeans(gIX_seed,M);
    % save motor seed functional trace, e.g. new motor regressors
    M_motorseedRegs{i_fish} = Reg;
    
    %% save motor seeds in VAR, and visualize (not in tRes)
    clusgroupID = 11;
    clusIDoverride = 2;
    name = 'Motor_seed_030317';
    SaveCluster_Direct(cIX_seed,gIX_seed,absIX,i_fish,name,clusgroupID,clusIDoverride);

    %% left-right combined plot
    figure('Position',[50,100,1400,800]);
    % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
    subplot(121)
    setappdata(hfig,'isPlotBehavior',1);
    setappdata(hfig,'isStimAvr',0);
    setappdata(hfig,'isPlotLines',0);
    UpdateTimeIndex(hfig);
    DrawTimeSeries(hfig,cIX_seed,gIX_seed);
    
    % right plot
    ax = subplot(122)
    I = LoadCurrentFishForAnatPlot(hfig,cIX_seed,gIX_seed);
    DrawCellsOnAnat(I,ax);
    
    if isSaveFig,
        filename = fullfile(saveDir2, ['Fish',num2str(i_fish),'_motormap_seed']);
        saveas(gcf, filename, 'png');
        close(gcf)
    end

    %% regression (in tRes), thresholding by % of cells (instead of corr thres)
    
    [~,Reg_tRes] = GetTrialAvrLongTrace(hfig,Reg);
    
    Corr = corr(Reg_tRes',M_0_tRes');
    [corr_max,IX] = max(Corr,[],1);
    [~,I] = sort(corr_max,'descend');
    
    M_target_prctcell= [1,2,3];
    for i_numcell = 2,
        nCells_total = size(M_0,1);
        prctcell = M_target_prctcell(i_numcell);
        nCells_target = round(prctcell/100 * nCells_total);
        
        %% target cell number: counting cells in hindbrain only
        
        Msk_IDs = 114; % mask for full hindbrain
        
        % isScreenMskFromAllCells
        cIX = (1:length(absIX))';
        gIX = ones(size(cIX));
        [cIX_hb,gIX_hb] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);

        I_hb = ismember(I,cIX_hb);
        cum_I_hb = cumsum(I_hb);
        lastIX = find(cum_I_hb==nCells_target);
        
        % cut off the desired number of cells to display
        cIX = I(1:lastIX)';
        gIX = IX(cIX)';
        M = UpdateIndices_Manual(hfig,cIX,gIX);
        
        M_thres_reg(i_numcell,i_fish) = corr_max(I(nCells_target));
        M_compareMotorCellNumber(1,i_fish) = length(cIX);
        
        %% k-means (not in tRes)
        numK = 16; % manual
        gIX = Kmeans_Direct(M,numK);
        
        %% Plot anat
        if isPlotFig,
            %             figure('Position',[50,100,800,1000]);
            %             I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
            %             DrawCellsOnAnat(I);
            figure('Position',[50,100,1400,800]);
            % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
            subplot(121)
            setappdata(hfig,'isPlotBehavior',1);
            setappdata(hfig,'isStimAvr',0);
            setappdata(hfig,'isPlotLines',0);
            UpdateTimeIndex(hfig);
            DrawTimeSeries(hfig,cIX,gIX);
            
            % right plot
            ax = subplot(122)
            I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
            DrawCellsOnAnat(I,ax);

            %% Save plot
            if isSaveFig,
                filename = fullfile(saveDir1, ['Fish',num2str(i_fish),'_motormap_',num2str(prctcell),'%']);
                saveas(gcf, filename, 'png');
                close(gcf)
            end
        end
        %%
%         %% Rank by stim-lock
%         [gIX,rankscore] = RankByStimLock_Direct(hfig,gIX);
%         
%         setappdata(hfig,'rankscore',round(rankscore*100)/100);
%         setappdata(hfig,'clrmap_name','jet');
%         
%         % Plot anat
%         if isPlotFig,
%             figure('Position',[50,100,800,1000]);
%             I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
%             DrawCellsOnAnat(I);
%             
%             
%             % Save plot
%             if isSaveFig,
%                 filename = fullfile(saveDir2, ['Fish',num2str(i_fish),'_motormap_',num2str(target_numcell),'_stimlock_anat']);
%                 saveas(gcf, filename, 'png');
%                 close(gcf)
%             end
%             
%             % Plot functional traces
%             figure('Position',[50,100,800,1000]);
%             % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
%             setappdata(hfig,'isPlotBehavior',1);
%             setappdata(hfig,'isStimAvr',0);
%             UpdateTimeIndex(hfig);
%             DrawTimeSeries(hfig,cIX,gIX);
%             
%             % Save plot
%             isSaveFig
%             filename = fullfile(saveDir2, ['Fish',num2str(i_fish),'_motormap_',num2str(target_numcell),'_stimlock_func']);
%             saveas(gcf, filename, 'png');
%             close(gcf)
%         end
%     end
%     % reset colormap
%     setappdata(hfig,'clrmap_name','hsv_new');


    end
    
    %% count cell numbers regressed (in tRes) with fictive trace
    thres_reg = M_thres_reg(2,i_fish);
%     [~,~,regressor_m_raw] = GetMotorRegressor(behavior,i_fish);
    
    [rawmotorcorr,IX] = max(corr(regressor_m_raw_tRes',M_0_tRes'),[],1);
    ix = find(rawmotorcorr>thres_reg);
    M_compareMotorCellNumber(2,i_fish) = length(ix);
end
SaveVARwithBackup();

%%
% range_fish excludes Fish 4
M_compareMotorCellNumber(:,4) = NaN;
figure;bar(M_compareMotorCellNumber')
