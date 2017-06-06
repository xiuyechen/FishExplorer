% fig2: anat averages from paired regressors: sensory/motor correlations
% plots a pair of regressors in a graded double colormap, saves as tiff
% stacks, and also saves anat average plot as tiff file.
% (compare to fig2A_correlation_1reg_with_hist for more details)

clear all; close all; clc

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% set params
thres_reg_const = 0.5;
% M_reg_thres = {0.5,0.5,0.5};

isKeepTopPrct = 0; % init


for setflag = 8
    % reset
    ResetDisplayParams(hfig);
    
    switch setflag
        case 1 % Phototaxis
%             stimrange = 1;
M_stimrange = GetStimRange('P');
            M_regtypeflag = [1,1,2]; % 0 for VAR, 1 for stim and 2 for motor
            M_reg_name = {'PT','HalfField','PT_SwimLR'};
            M_reg_range = {[3,2],[5,6],[1,3]};
            M_fishrange = {[1:18],[1:7],[1:3,5:18]};
            n_reg = length(M_reg_name);
            
        case 2 % OMR
%             stimrange = 2;
M_stimrange = GetStimRange('O');
            M_regtypeflag = [1,1,2]; % 1 for stim and 0 for motor
            M_reg_name = {'OMR_FwBw','OMR_LR','OMR_SwimLR'};
            M_reg_range =  {[7,6],[9,8],[1,3]};
            M_fishrange = {[8:18],[8:18],[8:18]};
            n_reg = length(M_reg_name);
            
        case 3 % Looming
%             stimrange = 5;
M_stimrange = GetStimRange('L');
            M_regtypeflag = [1,2]; % 1 for stim and 0 for motor
            M_reg_name = {'Loom_LR','Loom_SwimLR'};
            M_reg_range =  {[11,12],[1,3]};
            M_fishrange = {[9:15,17:18],[9:15,17:18]};
            n_reg = length(M_reg_name);
            
        case 4 % Dark Flash (black-white)
%             stimrange = 3;
M_stimrange = GetStimRange('D');
            M_regtypeflag = [1,2]; % 1 for stim and 0 for motor
            M_reg_name = {'DF_BW','DF_SwimLR'};
            M_reg_range =  {[1,4],[1,3]};
            M_fishrange = {[12:15,17:18],[12:15,17:18]};
            n_reg = length(M_reg_name);
            
%             % from VAR:

        case 5 % ABN
            M_regtypeflag = [0];
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'ABN_reg0.5_defstimrange_test'};
            M_clus_range = {[12,1]};
            M_fishrange = {[1:12,14:18]};
            n_reg = length(M_reg_name);
        case 6 % Fw
            M_regtypeflag = [0];
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'Fw_reg0.5_defS'}; 
            M_clus_range = {[11,4]};
            M_fishrange = {[1:18]};
            n_reg = length(M_reg_name);
            
        case 7 % HBO 4 stripes
            M_regtypeflag = [0,0];
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'HBO-L2_reg1%_tRes_defS','HBO-R2_reg1%_tRes_defS'}; 
            M_clus_range = {[10,2],[10,3]};
            M_fishrange = {[1:3,5:18],[1:3,5:18]};
            n_reg = length(M_reg_name);
            isKeepTopPrct = 1;
            prct_const = 1;
            setappdata(hfig,'isTrialRes',1);
            
        case 8 % HBO 2 halves
            M_regtypeflag = [0,0];
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'HBO2_reg1%_tRes_defS'}; 
            M_clus_range = {[10,5]};
            M_fishrange = {[1:3,5:18]};
            n_reg = length(M_reg_name);
            isKeepTopPrct = 1;
            prct_const = 1;
            setappdata(hfig,'isTrialRes',1);
    end
    M_fishrange_im = M_fishrange;
    
    %% Load fish
    range = 1:18;
    IM_full = cell(n_reg,max(range));
    %%
    
%     [M_fishrange_im,fishrange_load] = CheckIfLoadFish(M_fishrange,M_ClusterIDs,M_stimrange);
% [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,M_ClusterIDs{i_set},M_stimrange{i_set});
    for i_fish = range
        if ismember(i_fish,cell2mat(M_fishrange))
            ClusterIDs = [1,1];%GetClusterIDs('all');
            stimrange = M_stimrange{i_fish};   
            % load fish data
            if isempty(stimrange)
                [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
            else
                [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
            end
        end
        
        %% Load stim/motor
        for i_set = 1:n_reg
            
            if ~ismember(i_fish,M_fishrange{i_set})
                continue;
            end
            
            reg_thres = thres_reg_const; %M_reg_thres{i_set};
            
            %% Get Regressors
            if M_regtypeflag(i_set)==0
                ClusterIDs = M_clus_range{i_set};%[12,1];% GetClusterIDs('all');
                [cIX,gIX] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2));
                M = UpdateIndices_Manual(hfig,cIX,gIX);
                Reg = FindClustermeans(gIX,M);
                
            elseif M_regtypeflag(i_set)==1
                fishset = getappdata(hfig,'fishset');
                [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
                        
                reg_range = M_reg_range{i_set}; % left/right pair
                Reg = regressors(reg_range,:);
                            
            elseif M_regtypeflag(i_set)==2
                isMotorseed = 0;
                setappdata(hfig,'isMotorseed',isMotorseed);
                [~,~,behavior] = UpdateTimeIndex(hfig);                
                [~,names,regressors] = GetMotorRegressor(behavior,i_fish);           
                            
                reg_range = M_reg_range{i_set}; % left/right pair
                Reg = regressors(reg_range,:);
            end
            
            %% Regression
            % code adapted from 'best regressor regression' code 'AllRegsRegression'

            Corr = corr(Reg',M_0');
            
            if isKeepTopPrct
                % keep best regression only
                [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
                
                % top x percentile
                nCells_total = size(M_0,1);
%                 prct_const = 1;
                [CIX,RegThres] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_rows);
                cIX1 = CIX{1};
                cIX2 = CIX{2};
                
                % get map color
                clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres,1],64);
                clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres,1],64);
                %     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres1,1],64);
                %     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres2,1],64);
                
                cIX = [cIX1;cIX2];
                %     if isempty(cIX)
                %         M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish); %#ok<SAGROW>
                %         continue;
                %     end
                
                clrIX = [clrIX1;clrIX2];
                gIX_offset = [ones(size(cIX1));2*ones(size(cIX2))];
                gIX = clrIX+(gIX_offset-1)*64;
                %     gIX = [clrIX1;64+clrIX2];
                numK = length(unique(gIX));
                
            else
                %%
                [corr_max,IX_regtype] = max(Corr,[],1);
                cIX = find(corr_max>reg_thres)';
                gIX_offset = IX_regtype(cIX)';
                clrIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
                
                if isempty(cIX)
                    M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish);
                    continue;
                end
                
                gIX = clrIX+(gIX_offset-1)*64;
                numK = length(unique(gIX));
            end
            
            %% make double colormap
            clr1 = [1,0,0];
            clr1_ = [0.7,0.5,0.5];
            clr2 = [0,1,1];
            clr2_ = [0.5,0.7,0.7];
            numC = 64;
            clrmap1 = Make1DColormap([clr1_;clr1],numC);
            clrmap2 = Make1DColormap([clr2_;clr2],numC);
            clrmap = [clrmap1;clrmap2];
            
            %% make figure
            I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
            [h,im_full] = DrawCellsOnAnat(I);
            
            %         % add 2 colorbars
            % %         AddColorbarToAnat(clrmap,cmin,cmax)
            % %         colormap(clrmap1);
            % %         caxis([reg_thres,1])
            % %         colorbar('Location','manual','Position',[0.8,0.7,0.05,0.15],'Units','normalized')
            %         ax = axes('Position',[0.75,0.8,0.05,0.15],'Units','normalized');
            %         DrawCustomColorbar(clrmap1,[reg_thres,1],2,ax);
            %
            %         ax = axes('Position',[0.9,0.8,0.05,0.15],'Units','normalized');
            %         DrawCustomColorbar(clrmap2,[reg_thres,1],2,ax);
            %
            %% save figure
            close(h);
            IM_full{i_set,i_fish} = im_full;
            
        end
    end
    
    %% save as tiff stack
    for i_set = 1:n_reg
        range_im = M_fishrange_im{i_set};
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
        IM = IM_full(i_set,range_im);
        
        SaveImToTiffStack(IM,tiffdir);
    end
    
    %% [for later] plot from tiff stack
    isPlotfromtiffstack = 0;
    
    if isPlotfromtiffstack
        IM_full = cell(n_reg,18);
        for i_set = 1:n_reg
            %% get tiff-stack path
            tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
            
            %% load
            
            for i_fish = 1:18
                im = double(imread(tiffdir,i_fish))./255;
                IM_full{i_set,i_fish} = im(317:1236,1:621,:);
            end
        end
    end
    
    %% make average plots
    % M_k_scale = {1,1.5,1};
    % M_k_contrast = {1.2,1.5,1.2};
    
    for i_set = 1:n_reg
        range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
        cellarray = IM_full(i_set,range_im);
        
        % adjust params for visualization
        k_scale = 0.4;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
end
