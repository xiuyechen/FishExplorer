% 2D plot, PT vs OMR stim corr: multistim, similar to best stim reg



% this script performs two sets of regression (second set with both tRes
% and lrRes), and calculates 'setdiff' in both directions, respectively.
% i.e. A minus B, B minus A. These two anat plots are then concatenated.

clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;


ClusterIDs = [11,2]; % init; can overrride
prct_const = 2; % init; can overrride
range_fish = 8:18;%18;%GetFishRange;% init; can overrride

caseflag = 2;
switch caseflag % NOTE: regressors hard-coded!
    %     case 1
    %         isSetDiffnotIntersect = 1;
    %         is1RegInA = 0;
    %         M_isTrialRes = [0,1];
    %         M_reg_name{1} = 'motormap_tAvr_VS_tRes_motorseed2_setdiff';
    
    case 1 % for fig6a
        isSetDiffnotIntersect = 0;
        is1RegInA = 0;
        M_isTrialRes = [0,0];
        M_reg_name{1} = 'PTintOMR_regbased';%'stimmaps_PT_VS_OMR_intunion';
        M_reg_range = {[3,2],[9,8]};
          
    case 2 % for fig6a
        isSetDiffnotIntersect = 0;
        is1RegInA = 0;
        M_isTrialRes = [0,0];
        M_reg_name{1} = 'OMRintLm_regbased';%'stimmaps_PT_VS_OMR_intunion';
        M_reg_range = {[9,8],[11,12]};
        range_fish = [9:15,17:18];

        %     case 3 % not used
        %         isSetDiffnotIntersect = 1;
        %         is1RegInA = 1;
        %         M_isTrialRes = [1,1];
        %         M_reg_name{1} = 'motormap_lrAND_VS_lrtRes_motorseed2_setdiff';
        % %         M_reg_name{1} = 'motormap_lrtAvr_VS_lrtRes_motorseed2_setdiff';
        %
        %     case 4 % good for forward seed
        %         isSetDiffnotIntersect = 1;
        %         is1RegInA = 1;
        %         M_isTrialRes = [0,0];
        %         M_reg_name{1} = 'motormap_lrAvr_VS_lrRes_motorseed2_setdiff';
        %
        %     case 5 % good for AHC/HBO flip? setdiff prob not worth it
        %         isSetDiffnotIntersect = 1;
        %         is1RegInA = 0;
        %         M_isTrialRes = [0,1];
        %         M_reg_name{1} = 'motormap_df_VS_lrtRes_motorseed2_setdiff';
        %
        %     case 6 % eye map! main
        %         isSetDiffnotIntersect = 1;
        %         is1RegInA = 0;
        %         M_isTrialRes = [0,1];
        %         M_reg_name{1} = 'eyemap_tAvr_VS_tRes_setdiff';
        %
        %         ClusterIDs = [12,1]; % override
        %         prct_const = 0.5; % override
        %         range_fish = [1:8,11,12,14:17]; % skip 13, difficult 9,10,18
end

% stimrange = 2;

M_fishrange_im{1} = range_fish;
n_reg = 1;

i_set = 1;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',1);


%% run fish
range_prct = 1:10;
Intersect_cIX = cell(18,length(range_prct));
Intersect_gIX = cell(18,length(range_prct));

M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

IM_full = cell(n_reg,18);
IM_AB = cell(n_reg,18);

M_hb_count = zeros(18,6);

%%
for i_fish = range_fish
    
    
    %% regression (in tRes), thresholding by % of cells (instead of corr thres)
    %     M_cIX = cell(1,2);
    %     M_gIX = cell(1,2);
    M_cIX = cell(length(range_prct),2);
    M_gIX = cell(length(range_prct),2);
    for i_itr = 1:2
        setappdata(hfig,'isTrialRes',M_isTrialRes(i_itr));
        
        [cIX_seed,gIX_seed,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
        
        %         Reg = FindClustermeans(gIX_seed,M);
        fishset = getappdata(hfig,'fishset');
        [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
        
        reg_range = M_reg_range{i_itr}; % left/right pair
        Reg = regressors(reg_range,:);
        
        %% compute correlation
        Corr = corr(Reg',M_0');
        
        % keep best regression only
        [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
        
        for i_prct_count = 1:length(range_prct)
            prct_const = range_prct(i_prct_count);
            % top % (used 2% for original fig6a)
            nCells_total = size(M_0,1);
            [CIX,RegThres] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_rows);
            
            if i_itr == 1 && is1RegInA
                cIX = CIX{1};
                % get map color
                reg_thres = 0.25;
                gIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
            else
                %%
                cIX1 = CIX{1};
                cIX2 = CIX{2};
                
                % get map color
                reg_thres = 0.25;
                clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres,1],64);
                clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres,1],64);
                %     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres1,1],64);
                %     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres2,1],64);
                
                cIX = [cIX1;cIX2];
                
                % get gIX
                clrIX = [clrIX1;clrIX2];
                gIX_offset = [ones(size(cIX1));2*ones(size(cIX2))];
                gIX = clrIX+(gIX_offset-1)*64;
                %     gIX = [clrIX1;64+clrIX2];
                numK = length(unique(gIX));
                
            end

            M_cIX{i_prct_count,i_itr} = cIX;
            M_gIX{i_prct_count,i_itr} = gIX;
        end
    end % i_itr ~ comparison

    %% intersection
    for i_prct_count = 1:length(range_prct)
        prct_const = range_prct(i_prct_count);
        
        % get cIX
        cIX1 = M_cIX{i_prct_count,1};
        cIX2 = M_cIX{i_prct_count,2};
        [cIX_int,ix] = intersect(cIX1,cIX2);        
        Intersect_cIX{i_fish,i_prct_count} = cIX_int;
        
        % get gIX
        gIX1 = M_gIX{i_prct_count,1};
        gIX_int = gIX(ix);       
        Intersect_gIX{i_fish,i_prct_count} = gIX_int;
    end
end

save(fullfile(outputDir,[M_reg_name{1},'_sweepthres.mat']),'Intersect_cIX','Intersect_gIX');
