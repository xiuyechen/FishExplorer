% 2D plot, PT vs OMR stim corr: multistim, similar to best stim reg


clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));


ClusterIDs = [2,1];%[11,2]; % init; can overrride
prct_const = 1;% til 12/22/17; % init; can overrride
range_fish = 8:18;%18;%GetFishRange;% init; can overrride

caseflag = 1;
switch caseflag % NOTE: regressors hard-coded!
    case 1 % downstream of fig6a
        isSetDiffnotIntersect = 0;
        is1RegInA = 0;
        M_isTrialRes = [0,0];
        M_reg_name{1} = 'PTintOMR_MO_period-thres';%'PTintOMR_MO_SMTthres';
        M_reg_range = {[3,2],[9,8]};
        M_stimrange = {1,2};
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
IM_2 = cell(2,18); % anat map left/right (w convergence cells)
IM_3 = cell(2,18); % anat map left/right
IM_int = cell(1,18); % for intersection

PTintOMR = cell(2,18);% left pair, right pair

M_hb_count = zeros(18,6);

%%
for i_fish = range_fish
    
    
    %% regression (in tRes), thresholding by % of cells (instead of corr thres)
    M_cIX = cell(1,2);
    M_gIX = cell(1,2);
    for i_itr = 1:2
        setappdata(hfig,'isTrialRes',M_isTrialRes(i_itr));
        
        [cIX_all,gIX_all,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,M_stimrange{i_itr});
        
        %         Reg = FindClustermeans(gIX_seed,M);
        fishset = getappdata(hfig,'fishset');
        [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
        
        reg_range = M_reg_range{i_itr}; % left/right pair
        Reg = regressors(reg_range,:);
        
        %% compute correlation
        Corr = corr(Reg',M_0');
        
        % keep best regression only
        [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
        
        % top 2 %
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
            
            clrIX = [clrIX1;clrIX2];
            gIX_offset = [ones(size(cIX1));2*ones(size(cIX2))];
            gIX = clrIX+(gIX_offset-1)*64;
            %     gIX = [clrIX1;64+clrIX2];
            numK = length(unique(gIX));
        end
        %% pool stats
        %     M_thres_reg(1,i_fish) = reg_thres1;
        %     M_thres_reg(2,i_fish) = reg_thres2;
        %     M_compareMotorCellNumber(1,i_fish) = length(cIX1);
        %     M_compareMotorCellNumber(2,i_fish) = length(cIX2);
        M_cIX{i_itr} = cIX;
        M_gIX{i_itr} = gIX;
        
    end % i_itr ~ comparison
    
    %% Section 1: make the setdiff/intersection plots
    
    %% make double colormap (for intersection map)
    clr1 = [1,0,0];
    clr1_ = [0.5,0.4,0.4];
    %     clr1_ = [0.7,0.5,0.5];
    clr2 = [0,1,1];
    clr2_ = [0.4,0.5,0.5];
    %     clr2_ = [0.5,0.7,0.7];
    numC = 64;
    clrmap1 = Make1DColormap([clr1_;clr1],numC);
    clrmap2 = Make1DColormap([clr2_;clr2],numC);
    clrmap = [clrmap1;clrmap2];
    
    %% intersection: PT & OMR
    [cIX,ix] = intersect(M_cIX{1},M_cIX{2});
    gIX = M_gIX{1}(ix);
    %     [set2,ix2] = setdiff(CIX{2},CIX{1});
    %     cIX = [set1;set2];
    %     gIX = [GIX{1}(ix1);GIX{2}(ix2)];
    % make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,~,im] = DrawCellsOnAnat(I);
    close(h);
    IM_int{i_fish} = im;
    
    % save cells
    c1 = SelectClusterRange(cIX,gIX,1:64);
    c2 = SelectClusterRange(cIX,gIX,65:128);
    PTintOMR{1,i_fish} = c1;
    PTintOMR{2,i_fish} = c2;
    
    %% main loop for left/right motor
    for i_lr = 1:2
        %% scatter plot with 4D components
        % get betas for this fish
        
        betas = Betas{i_lr,i_fish};
        b1 = betas(:,1);
        b2 = betas(:,2);
        b3 = betas(:,3);
        
        % set up plot dimensions
        caseflag = 2;
        switch caseflag
            case 1
                X = b1;
                Y = b2;
                Xname = 'motor only (b1)';
                Yname = 'SMT (b2)';
            case 2
                X = b1;%b3;%b1;
                Y = sqrt(b2.^2+b3.^2);%b2;
                Xname = 'motor only (b1)';
                Yname = 'periodic';
        end
        
        numcell = length(b1);
        
        A = X;
        topN = round(0.01*numcell); % top 5% cutoff
        [~,IX] = sort(A,'descend');
        thresA = A(IX(topN));
        
        B = Y;
        topN = round(0.01*numcell); % top 5% cutoff
        [~,IX] = sort(B,'descend');
        thresB = B(IX(topN));

        IX_pass = union(find(A>=thresA),find(B>=thresB));
        IX_fail = intersect(find(A<thresA),find(B<thresB));%find(A<thresA);
                
        % get min/max
        x0 = min(X(IX_pass));
        x1 = max(X(IX_pass));
        y0 = min(Y(IX_pass));
        y1 = max(Y(IX_pass));
        
        gIX_in = (1:length(X))';
        
        %% make custom 2-D colormap
        %     grid = Make4color2Dcolormap;
        if true
            grid = MakeDiagonal2Dcolormap;
        end
        
        % map data to 2D colormap, and cluster sizes
        clrmap0 = MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
        clrmap = clrmap0;
        clrmap(IX_fail,:) = ones(length(IX_fail),3)*0.5;
        
        if length(clrmap)>1000
            U_size = ones(size(X));
        else % for actual clusters
            % get cluster size (number of cells in each cluster)
            gIX2 = SqueezeGroupIX(gIX_in);
            U = unique(gIX2);
            U_size = zeros(size(X));
            for i = 1:length(U)
                ix = find(gIX2 == U(i));
                U_size(i) = length(ix);
            end
        end
        
        %% bubble plot
        h = figure('Position',[500,500,300,250]); hold on
        scatter(X,Y,U_size,clrmap,'filled')
%         plot([x0,x1],[y0,y0],'k--');
%         plot([x0,x0],[y0,y1],'k--');
%         plot([x0,x1],[y0,y0],'k--');
        
        xlabel(Xname);ylabel(Yname);
        axis equal
        xlim([-0.5,1]);
        ylim([-0.5,1]);
        
        % highlight intersection cells
        [c1,g1] = SelectClusterRange(cIX,gIX,1:64);
        [c2,g2] = SelectClusterRange(cIX,gIX,65:128);
        if i_lr ==1
            c = c1;
            g = g1;
            clr = [1,0,0];
        else
            c = c2;
            g = g2;
            clr = [0,1,1];
        end
        
        if i_lr ==1
            scatter(X(c1),Y(c1),4,[1,0.5,0.5]);%,'filled');
        else
            scatter(X(c2),Y(c2),4,[0.5,1,1]);
        end
        
        % highlight anterior hindbrain (Rh1&2) cells (Rh1 219; Rh2 220;)
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        [c_hb,g_hb] = ScreenCellsWithMasks([219,220],c,g,MASKs,CellXYZ_norm,absIX);
        scatter(X(c_hb),Y(c_hb),4,clr,'filled');
        
        %% save bubble plot
        IM_1{i_lr,i_fish} = print('-RGBImage');
        close(h)
        
        %% anat map of bubble plot pass-y
        % highlight the intersection cells
        clrmap_plot = clrmap;
        clrmap_plot(g,:) = repmat(clr,length(g),1);
        
        cIX_plot = [cIX_all(IX_pass);c];
        gIX_plot = [gIX_in(IX_pass);g];
        I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap_plot);
        [h,~,im] = DrawCellsOnAnat(I);
        
        %% save
        close(h)
        IM_2{i_lr,i_fish} = im;
        
        %% just SM anat map
        % don't highlight the intersection cells
        clrmap_plot = clrmap;
        
        cIX_plot = [cIX_all(IX_pass)];
        gIX_plot = [gIX_in(IX_pass)];
        I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap_plot);
        [h,~,im] = DrawCellsOnAnat(I);
        
        %% save
        close(h)
        IM_3{i_lr,i_fish} = im;
    end
    
end


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
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_anat_allfish.tiff']);
        IM = IM_2(i_lr,range_im);
        SaveImToTiffStack(IM,tiffdir);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_anat_wo-int_allfish.tiff']);
        IM = IM_3(i_lr,range_im);
        SaveImToTiffStack(IM,tiffdir);
    end
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_allfish.tiff']);
    IM = IM_int(range_im);
    SaveImToTiffStack(IM,tiffdir);
end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

for i_set = 1:n_reg
    range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    %    L: [8,10,12,13,17]
    %    R: [8,9,11,16]
    for i_lr = 1:2
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
        cellarray = IM_2(i_lr,range_im);
        
        % adjust params for visualization
        k_scale = 0.6;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1.5;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_anat_avr_select5.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
    %% intersection map (uses individual stimrange, while fig6a used full range
    cellarray = IM_int(range_im);
    
    % adjust params for visualization
    k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.2;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_int_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end
