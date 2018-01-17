% 12/26/17: This code was inherited and developed from several setdiff/4D
% scripts. The 4D decomposition idea applied to find multi-stim convergence
% cells. Plot # of convergence cells as a function of % cutoff.

clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;

ClusterIDs = [2,1];%[11,2]; % init; can overrride
% prct_const = 3;% til 12/22/17; % init; can overrride


caseflag = 1;
switch caseflag % NOTE: regressors hard-coded!
    case 1 % downstream of fig6a, PT vs OMR
        load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));
        M_reg_name{1} = 'PTintOMR_MO_period-thres';
%         M_reg_range = {[3,2],[9,8]};
        M_stimrange = {1,2};
        range_fish = 8:18;
        range_prct = 1:10;
%         prct_const = 3;
        
    case 2 % OMR vs looming
        load(fullfile(outputDir,'4D_SM_stimrangeOMRloom_betas.mat'));
        M_reg_name{1} = 'OMRintLoom_MO_period-thres';
        M_stimrange = {2,5};
        range_fish = [9:15,17:18];
        prct_const = 8;
        
    case 3 % PT vs looming
        load(fullfile(outputDir,'4D_SM_stimrangePTloom_betas.mat'));
        M_reg_name{1} = 'PTintLoom_MO_period-thres';
        M_stimrange = {1,5};
        range_fish = [9:15,17:18];
        prct_const = 8;%3;
             
    case 4 % PT vs DF
        load(fullfile(outputDir,'4D_SM_stimrangePTDF_betas.mat'));
        M_reg_name{1} = 'PTintDF_MO_period-thres';
        M_stimrange = {1,3};
        range_fish = [12:15,17:18];
        prct_const = 2;   
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

i_lr = 2;

i_fish = 8;
% for i_fish = range_fish
    %% get top cells from individual stimrange (to do intersection later)
    M_cIX = cell(length(range_prct),2);
    M_cIX_diff = cell(length(range_prct),2);
    M_thresA = zeros(length(range_prct),1);
    M_motor_b1 = zeros(length(range_prct),2);
%     M_gIX = cell(length(range_prct),2);
    for i_itr = 1:2
        %%
        [cIX_all,gIX_all,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,M_stimrange{i_itr});
        
        
        [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace(hfig,M_0);

        [motor_tAvr,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace(hfig,behavior);
        %%
        b_stim = sqrt(abs(var(Data_tAvr')./var(Data_p')));
        b1_L = corr(motor_tRes(1,:)',Data_p');
        b1_R = corr(motor_tRes(2,:)',Data_p');
% %         b2 = corr(motor_tAvr(i_lr,:)',Data_p');
% %         b3 = sqrt(abs(var(Data_tAvr')./var(Data_p') - b2.^2)); % var(Data') is all 1
        
        % top %
        
        for i_prct_count = 1:length(range_prct)
            prct_const = range_prct(i_prct_count);
            A = b_stim;
            numcell = size(M_0,1);
            topN = round(prct_const/100*numcell); % top _% cutoff
            [A_sorted,IX] = sort(A,'descend');
            thresA = A_sorted(topN);
            
            IX_pass = find(A>thresA);            
            
            M_thresA(i_prct_count,1) = thresA;
            b1_LR = max([b1_L;b1_R],[],1);
            M_motor_b1(i_prct_count,1) = mean(b1_LR(IX_pass));
            M_motor_b1(i_prct_count,2) = max(b1_LR(IX_pass));
            M_cIX{i_prct_count,i_itr} = cIX_all(IX_pass);
            
            if i_prct_count==1
                M_cIX_diff{i_prct_count,i_itr} = cIX_all(IX_pass);
            else
                M_cIX_diff{i_prct_count,i_itr} = setdiff(cIX_all(IX_pass),M_cIX_diff{i_prct_count-1,i_itr});
            end
        end
    end % i_itr ~ comparison
    
    %% Section 1: make the setdiff/intersection plots
    
    %% make double colormap (for intersection map)
    %     clr1 = [1,0,0];
    %     clr1_ = [0.5,0.4,0.4];
    %     %     clr1_ = [0.7,0.5,0.5];
    %     clr2 = [0,1,1];
    %     clr2_ = [0.4,0.5,0.5];
    %     %     clr2_ = [0.5,0.7,0.7];
    %     numC = 64;
    %     clrmap1 = Make1DColormap([clr1_;clr1],numC);
    %     clrmap2 = Make1DColormap([clr2_;clr2],numC);
    %     clrmap = [clrmap1;clrmap2];
    
    %% intersection: PT & OMR
    %     [cIX_int,ix] = intersect(M_cIX{1},M_cIX{2});
    

    N_int = zeros(length(range_prct),1);
    cIX_grad = [];
    gIX_grad = [];
    for i_prct_count = 3:-1:1%length(range_prct)
        %%
        cIX1 = M_cIX{i_prct_count,1};
        cIX2 = M_cIX{i_prct_count,2};
        
        [cIX_int,ix] = intersect(cIX1,cIX2);
        gIX_int = ones(size(cIX_int));% not used here but for other plotting
        
        % get # of convergence cells as a function of % cutoff
        N_int(i_prct_count) = length(cIX_int);
        
        % make gIX index for full range of prct values
        gIX_grad = [gIX_grad;ones(size(cIX_int)).*i_prct_count];
        cIX_grad = [cIX_grad;cIX_int];
    end        

    % make figure
    setappdata(hfig,'clrmap_name','jet');
    I = LoadCurrentFishForAnatPlot(hfig,cIX_grad,gIX_grad);%,clrmap);
    [h,~,im] = DrawCellsOnAnat(I);
    %%
    close(h);
    IM_int{i_fish} = im;
    
    % save cells
%     Intersect_cIX{i_fish} = cIX_int;