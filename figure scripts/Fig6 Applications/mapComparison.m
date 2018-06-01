% map comparison test code

% fig2: anat averages from paired regressors: sensory/motor correlations
% plots a pair of regressors in a graded double colormap, saves as tiff
% stacks, and also saves anat average plot as tiff file.
% (compare to fig2A_correlation_1reg_with_hist for more details)
% function mapComparison()
clear all; close all; clc

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% set params
reg_thres = 0.5; % now only used for colormap maybe?? not for cell selection
% M_reg_thres = {0.5,0.5,0.5};

isKeepTopPrct = 1; % init, can override
prct_const = 2;

%% init
P_len_setA = zeros(4,4,18);
P_len_setB = zeros(4,4,18);
P_len_setAB = zeros(4,4,18);
IM_split = cell(4,4,18);
IM_3clr = cell(4,4,18);

M_stim_range = cell(1,2);
M_reg_name = cell(1,2);
M_reg_range = cell(1,2);

%% main loop
for ii = 1:3
    [M_stim_range{1},M_reg_name{1},M_reg_range{1},M_fishrange1] = getStimSetParams(ii);
    
    for jj = (ii+1):4
        [M_stim_range{2},M_reg_name{2},M_reg_range{2},M_fishrange2] = getStimSetParams(jj);

        % reset
        ResetDisplayParams(hfig);
        
        
        
        n_reg = 1;%length(M_reg_name);
        
        range_fish = intersect(M_fishrange1{:},M_fishrange2{:});
        M_fishrange = {range_fish};
        M_fishrange_im = M_fishrange;
        
        
        %         range = 1:18; % init,M_fishrange below specifies actual range #oldcode
        IM_AorB = cell(n_reg,18);
        
        P = zeros(18,3); % pool numbers
        
        %%
        
        
        for i_fish = range_fish
            %             if ismember(i_fish,cell2mat(M_fishrange))
            ClusterIDs = [1,1];%GetClusterIDs('all');
            
            %% sets
            
            CIX_set = cell(1,2);
            GIX_set = cell(1,2);
            for i_set = 1:2
                % load fish data
                stimrange = M_stim_range{i_set}{i_fish};
                [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
                
                %% Get stim regressors
                fishset = getappdata(hfig,'fishset');
                [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
                
                reg_range = M_reg_range{i_set}; % left/right pair
                Reg = regressors(reg_range,:);
                
                %% Regression
                % code adapted from 'best regressor regression' code 'AllRegsRegression'
                
                Corr = corr(Reg',M_0');
                
                
                % keep best regression only
                [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
                
                % top x percentile
                nCells_total = size(M_0,1);
                %                 prct_const = 1;
                [CIX,RegThres] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_rows);
                
                %                     if length(CIX)==1
                %                         cIX = CIX{1};
                %                         % get map color
                %                         reg_thres = 0.25;
                %                         gIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
                %                     else
                cIX1 = CIX{1};
                cIX2 = CIX{2};
                
                % get map color
                clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres,1],64);
                clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres,1],64);
                cIX = [cIX1;cIX2];
                
                clrIX = [clrIX1;clrIX2];
                
                CIX_set{i_set} = cIX;
                GIX_set{i_set} = clrIX;
            end
            
            %% count and save
            
            [cIX_a,ia] = setdiff(CIX_set{1},CIX_set{2});
            gIX_a = GIX_set{1}(ia);%ones(length(ia),1);
            [cIX_b,ia] = setdiff(CIX_set{2},CIX_set{1});
            gIX_b = GIX_set{2}(ia);%ones(length(ia),1);
            [cIX_ab,ia,ib] = intersect(CIX_set{1},CIX_set{2});
            gIX_ab = round(mean([GIX_set{1}(ia),GIX_set{2}(ib)],2));%ones(length(ia),1);            
            
            % pool
            P_len_setA(ii,jj,i_fish) = length(CIX_set{1});
            P_len_setB(ii,jj,i_fish) = length(CIX_set{2});
            P_len_setAB(ii,jj,i_fish) = length(cIX_ab);
            
            %% make triple colormap
            % channel 1 red, channel 2 green, intersection yellow
            
            clr1 = [1,0,0];
            clr1_ = [0.7,0.5,0.5];
            clr2 = [0,1,0];
            clr2_ = [0.5,0.7,0.5];
            clr12 = [1,1,0];
            clr12_ = [0.7,0.7,0.5];
            %                         clr12 = [0,0,1];
            %                         clr12_ = [0.5,0.5,0.7];
            numC = 64;
            clrmap1 = Make1DColormap([clr1_;clr1],numC);
            clrmap2 = Make1DColormap([clr2_;clr2],numC);
            clrmap12 = Make1DColormap([clr12_;clr12],numC);
            clrmap = [clrmap1;clrmap2;clrmap12];
            
            %% plot non-overlapping cells
            cIX = [cIX_a;cIX_b];
            clrIX = [gIX_a;gIX_b];
            gIX_offset = [ones(size(cIX_a));2*ones(size(cIX_b))];
            gIX = clrIX+(gIX_offset-1)*64;
            
            % make figure
            I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
            [h,im_full1] = DrawCellsOnAnat(I);
            close(h);
            
            %% plot only overlapping cells
            cIX = [cIX_ab];
            clrIX = [gIX_ab];
            gIX_offset = [3*ones(size(cIX_ab))];
            gIX = clrIX+(gIX_offset-1)*64;
            
            % make figure
            I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
            [h,im_full2] = DrawCellsOnAnat(I);
            close(h);

            %% plot 3 clrs together
            cIX = [cIX_a;cIX_b;cIX_ab];
            clrIX = [gIX_a;gIX_b;gIX_ab];
            gIX_offset = [ones(size(cIX_a));2*ones(size(cIX_b));3*ones(size(cIX_ab))];
            gIX = clrIX+(gIX_offset-1)*64;
            
            % make figure
            I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
            [h,im_full3] = DrawCellsOnAnat(I);
            close(h);
              
            %% save figs
            border = ones(size(im_full1,1),20,3);
            IM_split{ii,jj,i_fish} = horzcat(im_full1,border,im_full2);
            IM_3clr{ii,jj,i_fish} = im_full3;
            
        end
        
        %% save as tiff stack
        set_name = [M_reg_name{1},'_',M_reg_name{2}];
        
        tiffdir = fullfile(outputDir,[set_name,'_split_allfish.tiff']);
        IM = IM_split(ii,jj,range_fish);        
        SaveImToTiffStack(IM,tiffdir);

        tiffdir = fullfile(outputDir,[set_name,'_3clr_allfish.tiff']);
        IM = IM_3clr(ii,jj,range_fish);        
        SaveImToTiffStack(IM,tiffdir);
        
        %% make average plots
        
        cellarray = IM_split(ii,jj,range_fish);        
        k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1.1;%M_k_contrast{i_set};
        [~,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        tiffdir = fullfile(outputDir,[set_name,'split_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');        
        
        cellarray = IM_3clr(ii,jj,range_fish);        
        k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1.1;%M_k_contrast{i_set};
        [~,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        tiffdir = fullfile(outputDir,[set_name,'3clr_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
        
    end
end

%%
% P_len_setA = zeros(4,4,18);
% P_len_setB = zeros(4,4,18);
% P_len_setAB = zeros(4,4,18);
save(fullfile(outputDir,'mapcomparison_numbers.mat'),'P_len_setA','P_len_setB','P_len_setAB');

%%
ratio = P_len_setAB./P_len_setA;
im = nanmean(ratio,3);
figure;imagesc(im)

ratio2 = P_len_setAB./P_len_setB;
im2 = nanmean(ratio2,3);
