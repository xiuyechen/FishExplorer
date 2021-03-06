% 12/26/17: This code was inherited and developed from several setdiff/4D
% scripts. The 4D decomposition idea applied to the multi-stim convergence
% population, with lots of trial&error/options. The script as is is not 
% cleaned up with options, only kept the last one, which is plotting MO on
% the x axis and Periodic (sqrt(SO^2+SMT^2)) on the y axis. Colormap is
% kept to the simplest - 3 solid colors only for x/y/convergence (no
% gradients, hard cutoff based on number of convergence cells, which is
% the intersection of top 3% PT and OMR, respectively).

% outputs: 2D-scatter plot, anat map of top-MO and top-stim, with or 
% without convergence cell shown in third color, anat map of convergence 
% cells alone, and stack averages. For functional traces, I have some code
% at the end of the 'sandbox_downstreamoffig6a.m' that cleans up the
% colormap gIX, and then I import it into the GUI and separate the groups
% into left/right sides before plotting avr functional traces. For
% top-stim, need to further divide group into enough clusters (e.g. PT
% specific, OMR specific, and not-shown background specific). Plot
% functional traces in small numbers with code (instead of GUI export) to
% save vector files instead of pictures. 


clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;

ClusterIDs = [2,1];%[11,2]; % init; can overrride
prct_const = 3;% til 12/22/17; % init; can overrride
i_perct = 3; % init; can overrride

caseflag = 1;
switch caseflag % NOTE: regressors hard-coded!
    case 1 % downstream of fig6a, PT vs OMR]
        load(fullfile(outputDir,'4D_SM_stimrangePTOMR_minmax_betas.mat'));
%         load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));% up till 1/20/18
        M_reg_name{1} = 'PTintOMR_MO_period-3%subthres';
%         M_reg_range = {[3,2],[9,8]};
        M_stimrange = {1,2};
        range_fish = 8:18;
        prct_const = 3;% used for period based regression if ifLoadfromfile=0
        ifLoadfromfile = 1;
        i_perct = 3;
        stimrange = [1,2]; % new? not used
        
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

M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

IM_1 = cell(2,18); % scatter plot left/right
IM_2 = cell(1,18); % anat map left/right (w convergence cells)
IM_3 = cell(1,18); % anat map left/right
IM_int = cell(1,18); % for intersection
IM_int_raw = cell(1,18); % for intersection
% Intersect_cIX = cell(1,18);

M_pool = cell(3,18);

%%
if ifLoadfromfile
    
    load(fullfile(outputDir,'PTintOMR_regbased_sweepthres.mat'),'Intersect_cIX');
%     load(fullfile(outputDir,'PTintOMR_sweepthres'),'Intersect_cIX');
end


for i_fish = range_fish
    if ~ifLoadfromfile
    %% get top cells from individual stimrange (to do intersection later)
    M_cIX = cell(1,2);
    M_gIX = cell(1,2);
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
        topN = round(prct_const/100*numcell); % top _% cutoff
        [A_sorted,IX] = sort(A,'descend');
        thresA = A_sorted(topN);
        
        IX_pass = find(A>thresA);

        M_cIX{i_itr} = cIX_all(IX_pass);

    end % i_itr ~ comparison
    
    [cIX_int,ix] = intersect(M_cIX{1},M_cIX{2});
    gIX_int = ones(size(cIX_int));%M_gIX{1}(ix);
    
    else        
        cIX_all = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange,0);
        cIX_int_raw = Intersect_cIX{i_fish,i_perct};
        gIX_int_raw = ones(size(cIX_int_raw));
    end
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

    % save cells
%     Intersect_cIX{i_fish} = cIX_int;
%     PTintOMR{i_fish} = cIX_int;
    
    %% main loop for left/right motor
    PassX_2 = cell(1,2);
    PassY_2 = cell(1,2);
    for i_lr = 1:2
        %% scatter plot with 4D components
        % loaded betas for this fish for the combo data (including both stim)        
        betas = Betas{i_lr,i_fish};
        % set up plot dimensions
        X = betas(:,1);%b3;%b1;
        Y = betas(:,5);
        %                 Y = sqrt(b2.^2+b3.^2);%b2;
        Xname = 'motor only (b1)';
        Yname = 'periodic';
        
        numcell = length(X);
        
        A = X;
        topN = length(cIX_int_raw);%length(M_cIX{2});%%round(0.01*numcell); % top 5% cutoff
        [~,IX] = sort(A,'descend');
        thresA = A(IX(topN));
        
        B = Y;
        topN = length(cIX_int_raw);%length(M_cIX{2});%round(0.01*numcell); % top 5% cutoff
        [~,IX] = sort(B,'descend');
        thresB = B(IX(topN));%min(Y(cIX_int));%

        IX_passX = setdiff(find(A>=thresA),find(B>=thresB));
        IX_passY = setdiff(find(B>=thresB),find(A>=thresA));%find(B>=thresB);%
        IX_pass = union(find(A>=thresA),find(B>=thresB));
        IX_fail = intersect(find(A<thresA),find(B<thresB));%find(A<thresA);
                
        % get min/max
        x0 = min(X(IX_pass));
        x1 = max(X(IX_pass));
        y0 = min(Y(IX_pass));
        y1 = max(Y(IX_pass));
        
        gIX_in = (1:length(X))';
        
        PassX_2{i_lr} = IX_passX;
        PassY_2{i_lr} = IX_passY;
        
        cIX_int = setdiff(intersect(cIX_int_raw,IX_passY),IX_passX);
        
        %% make colormap
    
        clr1 = [0.3,0.7,0.2];
        clr2 = [0.1,0.3,0.9];%[0.3,0.8,1];
        clr_fail = [0.5,0.5,0.5];
        clr_int_raw = [1,0.2,0.2];%[1,0.1,0];
        clr_int = [0.7,0.2,0.2];
        clrmap = ones(numcell,3);%MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
        clrmap(IX_passX,:) = clr1.*ones(length(IX_passX),3);
        clrmap(IX_passY,:) = clr2.*ones(length(IX_passY),3);
        clrmap(IX_fail,:) = clr_fail.*ones(length(IX_fail),3);
        clrmap(cIX_int_raw,:) = clr_int_raw.*ones(length(cIX_int_raw),3);
        clrmap(cIX_int,:) = clr_int.*ones(length(cIX_int),3);
        
        %% bubble plot
        h = figure('Position',[500,100,300,250]); hold on
        U_size = ones(size(X));
        scatter(X,Y,U_size,clrmap,'filled')
        plot([x0,x1],[thresB,thresB],'k--');
        plot([thresA,thresA],[y0,y1],'k--');
        %         plot([x0,x1],[y0,y0],'k--');
        
        scatter(X(IX_passX),Y(IX_passX),1,clr1);%,'filled');
        scatter(X(IX_passY),Y(IX_passY),1,clr2);%,'filled');
        scatter(X(cIX_int_raw),Y(cIX_int_raw),1,clr_int_raw);%[1,0.5,0.5]);%,'filled');
        
        xlabel(Xname);ylabel(Yname);
        axis equal
        
        xlim([-0.4,0.8]);
        ylim([0,1]);
        
        scatter(X(cIX_int),Y(cIX_int),1,clr_int);%[1,0.5,0.5]);%,'filled');
        
        
        % highlight anterior hindbrain (Rh1&2) cells (Rh1 219; Rh2 220;)
%         MASKs = getappdata(hfig,'MASKs');
%         CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
%         absIX = getappdata(hfig,'absIX');
%         [c_hb,g_hb] = ScreenCellsWithMasks([219,220],cIX_int_raw,gIX_int_raw,MASKs,CellXYZ_norm,absIX);
%         scatter(X(c_hb),Y(c_hb),1,'k');%clr_int);
        
        %% save bubble plot
        IM_1{i_lr,i_fish} = print('-RGBImage');
        close(h)
        
        %% intersection anat map
        setappdata(hfig,'clrmap_name','hsv_old');
        gIX_int = ones(size(cIX_int));
        I = LoadCurrentFishForAnatPlot(hfig,cIX_int,gIX_int,clr_int);
        [h,~,im] = DrawCellsOnAnat(I);
        close(h);
        IM_int{i_fish} = im;
        
        I = LoadCurrentFishForAnatPlot(hfig,cIX_int_raw,gIX_int_raw,clr_int_raw);%,clrmap);
        [h,~,im] = DrawCellsOnAnat(I);
        close(h);
        IM_int_raw{i_fish} = im;


    end
    
    %% set colors for anat map
    %     IX_passX2 = union(PassX_2{1},PassX_2{2});
    clr1_1 = [0.5,1,0.4]; % clr1 = [0.5,1,0.5];
    clr1_2 = [0.5,1,0.4]; %[0.4,1,0.5];
    clrmap(PassX_2{1},:) = clr1_1.*ones(length(PassX_2{1}),3);
    clrmap(PassX_2{2},:) = clr1_2.*ones(length(PassX_2{2}),3);
    
%     clr2_1 = [0.7,0.4,1]; % clr1 = [0.5,1,0.5];
    clr2_2 = [0.2,0.4,1];
%     clrmap(PassY_2{1},:) = clr2_1.*ones(length(PassY_2{1}),3);
    clrmap(PassY_2{2},:) = clr2_2.*ones(length(PassY_2{2}),3);
    
    IX_pass_2 = union(union(PassX_2{1},PassX_2{2}),union(PassY_2{1},PassY_2{2}));
    
        %% anat map of bubble plot pass-y
        % highlight the intersection cells
        clrmap_plot = clrmap;
        
        clrmap_plot(gIX_int,:) = repmat(clr_int,length(gIX_int),1);
        
        cIX_plot = [cIX_all(IX_pass_2);cIX_int];
        gIX_plot = [gIX_in(IX_pass_2);gIX_int];

        
        % draw hindbrain Rh1&2 mask
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        setappdata(hfig,'isShowMasks',1);
        setappdata(hfig,'Msk_IDs',[219,220,221]);
        setappdata(hfig,'isShowMskOutline',1);
        I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap_plot);
        [h,~,im] = DrawCellsOnAnat(I);
        
        %% save
        close(h)
        IM_2{i_fish} = im;
        
        %% just SM anat map
        % don't highlight the intersection cells
        clrmap_plot = clrmap;
        
        cIX_plot = [cIX_all(IX_pass_2)];
        gIX_plot = [gIX_in(IX_pass_2)];
        I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap_plot);
        [h,~,im] = DrawCellsOnAnat(I);
        
        %% save
        close(h)
        IM_3{i_fish} = im;

        Intersect_cIX_thres{i_fish,i_perct} = cIX_int;
    
end

%%
save(fullfile(outputDir,'PTintOMR_regbased_sweepthres.mat'),'Intersect_cIX_thres','-append');

%% draw color bars - to save???
% figure
% ax = axes('Position',[0.75,0.8,0.05,0.15],'Units','normalized');
% DrawCustomColorbar(clrmap1,[reg_thres,1],2,ax);
%
% ax = axes('Position',[0.9,0.8,0.05,0.15],'Units','normalized');
% DrawCustomColorbar(clrmap2,[reg_thres,1],2,ax);

%% save as tiff stack
M_lr = {'-L','-R'};
n_reg = 1;
for i_set = 1:n_reg
    range_im = M_fishrange_im{i_set};
    
    for i_lr = 1:2
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_scatter_allfish.tiff']);
        IM = IM_1(i_lr,range_im);
        SaveImToTiffStack(IM,tiffdir);
    end
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_anat_allfish.tiff']);
    IM = IM_2(range_im);
    SaveImToTiffStack(IM,tiffdir);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_anat_wo-int_allfish.tiff']);
    IM = IM_3(range_im);
    SaveImToTiffStack(IM,tiffdir);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_allfish.tiff']);
    IM = IM_int(range_im);
    SaveImToTiffStack(IM,tiffdir);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_raw_allfish.tiff']);
    IM = IM_int_raw(range_im);
    SaveImToTiffStack(IM,tiffdir);
end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

for i_set = 1:n_reg
    range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    %    L: [8,10,12,13,17]
    %    R: [8,9,11,16]
%     for i_lr = 1:2
        %         %% scatter plot
        %         cellarray = IM_1(i_lr,range_im);
        %
        %         % adjust params for visualization
        %         k_scale = 1;%1/1.5;%M_k_scale{i_set};
        %         k_contrast = 0.6;%M_k_contrast{i_set};
        %
        %         [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        %
        %         tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_scatter_avr.tiff']);
        %         imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
        
        %% anat map
        cellarray = IM_2(range_im);%IM_2(i_lr,range_im);
        
        % adjust params for visualization
        k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1.5;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_anat_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
        
        % without int cells
        cellarray = IM_3(range_im);%IM_2(i_lr,range_im);
        
        % adjust params for visualization
        k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1.5;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_anat_woint_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
        
%     end
    
    %% intersection map (uses individual stimrange, while fig6a used full range
    cellarray = IM_int(range_im);
    
    % adjust params for visualization
    k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.2;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
        
    % int_raw
    cellarray = IM_int_raw(range_im);
    
    % adjust params for visualization
    k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.2;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_raw_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end

%%
% save(fullfile(outputDir,'PTintOMR_3%each.mat'),'PTintOMR');
