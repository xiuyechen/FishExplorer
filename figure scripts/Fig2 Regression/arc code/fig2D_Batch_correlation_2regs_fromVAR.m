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
global VAR;

%% set params
reg_thres = 0.5;
% M_reg_thres = {0.5,0.5,0.5};

for stimflag = 2
    % reset
    ResetDisplayParams(hfig);
    
    switch stimflag
        case 1 % ABN
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'ABN_reg0.5_defstimrange'}; % one at a time for now
            M_reg_clus = {[12,1]};
            M_fishrange = {[1:12,14:18]};
            n_reg = length(M_reg_name);
        case 2 % Fw
            M_stimrange = GetStimRange();%('2');
            M_reg_name = {'Fw_reg0.5_defS'}; % one at a time for now
            M_reg_clus = {[11,4]};
            M_fishrange = {[1:18]};
            n_reg = length(M_reg_name);
    end
    M_fishrange_im = M_fishrange;
    
    %% Load fish
    range = 1:18;
    IM_full = cell(n_reg,max(range));
    %%
    i_set = 1;
    
    for i_fish = range
        if ismember(i_fish,cell2mat(M_fishrange))
            ClusterIDs = M_reg_clus{i_set};
            stimrange = M_stimrange{i_fish};   
            % load fish data
            if ~isempty(stimrange)              
                
                % check that VAR is not empty
                nClus = length(VAR(i_fish).ClusGroup{ClusterIDs(1)});
                if nClus>= ClusterIDs(2) % ~~loose approx ~~ exist record
                    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
                else
                    M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish);
                    continue;
                end
            end            
        end
        
        %% regression
%         for i_set = 1%:n_reg
            
            if ~ismember(i_fish,M_fishrange{i_set})
                continue;
            end
            
            % code adapted from 'best regressor regression' code 'AllRegsRegression'
            Reg = FindClustermeans(gIX_load,M);
            
            if stimflag==2
                Reg = mean(Reg,1);
            end
            
            Corr = corr(Reg',M_0');
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
            
            %% make double colormap - ??
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
            
%         end
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
    
    %%
    % M_k_scale = {1,1.5,1};
    % M_k_contrast = {1.2,1.5,1.2};
    
    for i_set = 1:n_reg;
        range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
        cellarray = IM_full(i_set,range_im);
        
        % adjust params for visualization
        k_scale = 1;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
end
