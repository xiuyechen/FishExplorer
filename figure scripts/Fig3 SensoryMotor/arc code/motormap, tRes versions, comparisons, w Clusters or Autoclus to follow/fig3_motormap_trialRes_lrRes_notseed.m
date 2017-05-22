% new procedure: generate motor seed/map in trialRes space
% find top 100 cells correlating to motor, screen with Rh4+Rh5(+Rh6??) masks

% corr sweep (0.5-0.7) to get top %2 cells, then kmeans, save; rank by stim-lock, save
clear all; close all; clc

%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
saveDir = [];
saveDir{1} = fullfile(outputDir,'motor_map_notseed_lrRes_1_kmeans_030417');
saveDir{2} = fullfile(outputDir,'motor_map_notseed_lrRes_2_autoclus_030417');
if ~exist(saveDir{1}, 'dir'), mkdir(saveDir{1}), end;
if ~exist(saveDir{2}, 'dir'), mkdir(saveDir{2}), end;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',0);
setappdata(hfig,'isTrialRes',1);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];
M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

for i_fish = range_fish
    ClusterIDs = [1,1];
    [~,~,~,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');

    %% get motor regressors
    [~,~,regressor_m_raw] = GetMotorRegressor(behavior,i_fish);
    Reg = regressor_m_raw([1,3],:);
    
    %% regression (in tRes), thresholding by % of cells (instead of corr thres)
    
    %     [~,Reg_tRes] = GetTrialAvrLongTrace(hfig,Reg);
    %     [~,M_0_tRes] = GetTrialAvrLongTrace(hfig,M_0);
   
%     Reg = FindClustermeans(gIX_seed,M);
    Reg_LR = Reg-repmat(mean(Reg),2,1);
    
    %     Corr = corr(Reg_tRes_LR',M_0_tRes');
    Corr = corr(Reg_LR',M_0');
    [Corr_sorted,IX_corr] = sort(Corr,2,'descend');
    
    %     [corr_max,IX] = max(Corr,[],1);
    %     [~,I] = sort(corr_max,'descend');
    
    M_target_prctcell= [1,2,3];
    for i_numcell = 1:2,
        nCells_total = size(M_0,1);
        prctcell = M_target_prctcell(i_numcell);
        nCells_target = round(prctcell/100 * nCells_total);
        
        %% target cell number: counting cells in hindbrain only
        
        Msk_IDs = 114; % mask for full hindbrain
        
        % isScreenMskFromAllCells
        cIX = (1:length(absIX))';
        gIX = ones(size(cIX));
        [cIX_hb,gIX_hb] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
        
        I_hb = ismember(IX_corr,cIX_hb);
        cum_I_hb1 = cumsum(I_hb(1,:));
        cum_I_hb2 = cumsum(I_hb(2,:));
        lastIX{1} = find(cum_I_hb1==nCells_target,1,'first');
        lastIX{2} = find(cum_I_hb2==nCells_target,1,'first');
        
        % cut off the desired number of cells to display
        for i_lr = 1:2
            cIX = IX_corr(i_lr,1:lastIX{i_lr})';
            gIX = ceil((1:length(cIX))'/(length(cIX)/20));
            M = UpdateIndices_Manual(hfig,cIX,gIX);
            [~,M_tRes] = GetTrialAvrLongTrace(hfig,M);
            setappdata(hfig,'M',M_tRes);
            
            M_thres_reg(i_lr,i_fish) = Corr_sorted(i_lr,lastIX{i_lr});
            M_compareMotorCellNumber(i_lr,i_fish) = length(cIX);
            
            %% clustering: trying Autoclus here!
            if i_numcell==1
                %%
                numK = 16; % manual
                gIX = Kmeans_Direct(M_tRes,numK);
            else
                cIX_reg = cIX;%(1:size(M_0,1))';
                [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg);
                UpdateIndices_Manual(hfig,cIX,gIX);
            end
            %% Plot anat
            if isPlotFig,
            
                I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
                f{i_lr} = DrawCellsOnAnat(I);
%                 f{i_lr} = figure('Position',[50,100,1400,800]);
%                 % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
%                 subplot(121)
%                 setappdata(hfig,'isPlotBehavior',1);
%                 setappdata(hfig,'isStimAvr',0);
%                 setappdata(hfig,'isPlotLines',0);
%                 %             UpdateTimeIndex(hfig);
%                 DrawTimeSeries(hfig,cIX,gIX);
%                 
%                 % right plot
%                 ax = subplot(122);
%                 I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
%                 DrawCellsOnAnat(I,ax);

            end
        end
        
        combineFiguresLR(f{1},f{2});
        
        % Save plot
        if isSaveFig,
            filename = fullfile(saveDir{i_numcell}, ['Fish',num2str(i_fish),'_motormap_',num2str(prctcell),'%_each']);
            saveas(gcf, filename, 'png');
            close(gcf)
        end                        
        
    end
end

%%
% range_fish excludes Fish 4
% M_compareMotorCellNumber(:,4) = NaN;
% figure;bar(M_compareMotorCellNumber')
