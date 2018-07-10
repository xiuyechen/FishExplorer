% this script selects the top x% cells for each stim chunk, and finds the
% intersecting cells. gIX assigned with 1,2 for congruent pairs (same
% direction, e.g. PT-left and OMR-left) and 3,4 for incongruent pairs.
% sweeps percentage (for each regressor) from 1% to 10% and save as .mat

% (script inherited from setdiff scripts originally used for motor maps)

%%
clear all; close all; clc

%% init
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',1);

%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;

ClusterIDs = [11,2]; % doesn't matter % init; can overrride
% prct_const = 2; % init; can overrride
range_prct = 1:10;

setappdata(hfig,'isStimAvr',1);

% remnants
isSetDiffnotIntersect = 0;
is1RegInA = 0;
M_isTrialRes = [0,0];

for caseflag = 1:6
    switch caseflag
        % excluding fish10 for OMR (cong/incong confound)% 6/29/18
        case 1 % PT & OMR
            M_reg_name{1} = 'PTintOMR_regbased';
            M_reg_range = {[3,2],[9,8]};
            M_stimrange = {1,2};
            range_fish = [8,9,11:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
                            
        case 2 % PT & Looming
            M_reg_name{1} = 'PTintLm_regbased';
            M_reg_range = {[3,2],[11,12]};
            M_stimrange = {1,5};
            range_fish = [9:15,17:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
            
        case 3 % PT & DF
            M_reg_name{1} = 'PTintDF_regbased';
            M_reg_range = {[3,2],[1,4]}; % [3,1],[2,1] "congruent", [3,4],[2,4] "incongruent" % 3=PT-L, 1=DF
            M_stimrange = {1,3};
            range_fish = [1:5,12:15,17:18];%[12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]}; % PT&black, i.e. R-off, L-off
            M_incong_pairs = {[1,2],[2,2]}; % PT&white i.e. L-on, R-on
                
        case 4 % OMR & Looming
            M_reg_name{1} = 'OMRintLm_regbased';
            M_reg_range = {[9,8],[11,12]};
            M_stimrange = {2,5};
            range_fish = [9,11:15,17:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
            
        case 5 % OMR & DF
            M_reg_name{1} = 'OMRintDF_regbased';
            M_reg_range = {[9,8],[1,4]};
            M_stimrange = {2,3};
            range_fish = [12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]};
            M_incong_pairs = {[1,2],[2,2]};
            
        case 6 % Looming & DF
            M_reg_name{1} = 'LmintDF_regbased';
            M_reg_range = {[11,12],[1,4]};
            M_stimrange = {5,3};
            range_fish = [12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]};
            M_incong_pairs = {[1,2],[2,2]};
            
    end
    
    %% run fish
    
    Intersect_cIX = cell(18,length(range_prct));
    Intersect_gIX = cell(18,length(range_prct));
    
    n_reg = 1;
    IM_full = cell(n_reg,18);
    IM_AB = cell(n_reg,18);
    
    %%
    for i_fish = range_fish
        
        M_cIX = cell(length(range_prct),2,2); % dim2: PT or OMR; dim3: L or R (for DF: black or white)
        %     M_gIX = cell(length(range_prct),2,2);
        
        for i_itr = 1:2
            setappdata(hfig,'isTrialRes',M_isTrialRes(i_itr));
            
            if i_fish<8
                [cIX_seed,gIX_seed,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,1);
            else
                [cIX_seed,gIX_seed,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,M_stimrange{i_itr});
            end
            
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
                    %                 gIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
                else
                    %%
                    for i_reg = 1:2
                        cIX = CIX{i_reg};
                        
                        %                     % get map color
                        %                     reg_thres = 0.25;
                        %                     gIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
                        
                        M_cIX{i_prct_count,i_itr,i_reg} = cIX;
                        %                     M_gIX{i_prct_count,i_itr,i_reg} = gIX;
                        
                    end
                    
                end
            end
        end % i_itr ~ comparison
        
        %% intersection
        for i_prct_count = 1:length(range_prct)
            %         prct_const = range_prct(i_prct_count);
            
            % congruence/incongruence
            M_pairs = [M_cong_pairs,M_incong_pairs];
            
            cIX_int = [];
            gIX_int = [];
            for i_pairs = 1:length(M_pairs)
                pair = M_pairs{i_pairs};
                
                % get cIX
                cIX1 = M_cIX{i_prct_count,1,pair(1)};
                cIX2 = M_cIX{i_prct_count,2,pair(2)};
                [cIX_int0,ix] = intersect(cIX1,cIX2);
                gIX_int0 = ones(length(cIX_int0),1);
                
                cIX_int = [cIX_int;cIX_int0];
                gIX_int = [gIX_int;i_pairs*ones(length(cIX_int0),1)];
            end
            
            Intersect_cIX{i_fish,i_prct_count} = cIX_int;
            Intersect_gIX{i_fish,i_prct_count} = gIX_int;
        end
    end
    
    save(fullfile(outputDir,[M_reg_name{1},'_sweepthres_cong.mat']),'Intersect_cIX','Intersect_gIX');
    % save(fullfile(outputDir,[M_reg_name{1},'_sweepthres.mat']),'Intersect_cIX','Intersect_gIX');
end