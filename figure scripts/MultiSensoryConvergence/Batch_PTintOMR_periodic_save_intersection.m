
clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;

ClusterIDs = [2,1];%[11,2]; % init; can overrride
% prct_const = 3;% til 12/22/17; % init; can overrride
range_prct = 1:10;

caseflag = 2;
switch caseflag % NOTE: regressors hard-coded!
    case 1 % downstream of fig6a, PT vs OMR
        load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));
        M_reg_name{1} = 'PTintOMR';
        %         M_reg_range = {[3,2],[9,8]};
        M_stimrange = {1,2};
        range_fish = 8:18;
        
        
    case 2 % OMR vs looming
        load(fullfile(outputDir,'4D_SM_stimrangeOMRloom_betas.mat'));
        M_reg_name{1} = 'OMRintLoom_MO_period-thres';
        M_stimrange = {2,5};
        range_fish = [9:15,17:18];
        prct_const = 8;
        
        %     case 3 % PT vs looming
        %         load(fullfile(outputDir,'4D_SM_stimrangePTloom_betas.mat'));
        %         M_reg_name{1} = 'PTintLoom_MO_period-thres';
        %         M_stimrange = {1,5};
        %         range_fish = [9:15,17:18];
        %         prct_const = 8;%3;
        %
        %     case 4 % PT vs DF
        %         load(fullfile(outputDir,'4D_SM_stimrangePTDF_betas.mat'));
        %         M_reg_name{1} = 'PTintDF_MO_period-thres';
        %         M_stimrange = {1,3};
        %         range_fish = [12:15,17:18];
        %         prct_const = 2;
end

% stimrange = 2;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',1);


%% run fish

M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

IM_1 = cell(2,18); % scatter plot left/right
IM_2 = cell(1,18); % anat map left/right (w convergence cells)
IM_3 = cell(1,18); % anat map left/right
IM_int = cell(1,18); % for intersection

Intersect_cIX = cell(18,length(range_prct));

M_hb_count = zeros(18,6);

M_pool = cell(3,18);

%%
for i_fish = range_fish
    %% get top cells from individual stimrange (to do intersection later)
    M_cIX = cell(length(range_prct),2);
%     M_gIX = cell(length(range_prct),2);
    for i_itr = 1:2
        %%
        [cIX_all,gIX_all,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,M_stimrange{i_itr});
        
        
        [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace(hfig,M_0);
        
        [motor_tAvr,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace(hfig,behavior);
        
        b_stim = sqrt(abs(var(Data_tAvr')./var(Data_p')));
        % %         b1 = corr(motor_tRes(i_lr,:)',Data_p');
        % %         b2 = corr(motor_tAvr(i_lr,:)',Data_p');
        % %         b3 = sqrt(abs(var(Data_tAvr')./var(Data_p') - b2.^2)); % var(Data') is all 1
        
        % top %
        A = b_stim;
        numcell = size(M_0,1);
        
        for i_prct_count = 1:length(range_prct)
            prct_const = range_prct(i_prct_count);
            topN = round(prct_const/100*numcell); % top _% cutoff
            [A_sorted,IX] = sort(A,'descend');
            thresA = A_sorted(topN);
            
            IX_pass = find(A>thresA);
            
            M_cIX{i_prct_count,i_itr} = cIX_all(IX_pass);
        end
        
    end % i_itr ~ comparison
    
    %% intersection: PT & OMR    
    for i_prct_count = 1:length(range_prct)
        prct_const = range_prct(i_prct_count);
        
        cIX1 = M_cIX{i_prct_count,1};
        cIX2 = M_cIX{i_prct_count,2};
        
        [cIX_int,ix] = intersect(cIX1,cIX2);
        
        Intersect_cIX{i_fish,i_prct_count} = cIX_int;
    end
    
end

save(fullfile(outputDir,[M_reg_name{1},'_sweepthres.mat']),'Intersect_cIX');