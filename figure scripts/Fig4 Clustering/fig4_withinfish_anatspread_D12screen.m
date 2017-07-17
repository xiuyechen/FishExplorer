% heatmap: how much clusters in the same fish overlap anatomically
% distance metric same as previously used for between fish
clear all; close all; clc;

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range_fish = 1:18;
n_reg = 2;
IM_full = cell(n_reg,max(range_fish));

%%
% i_set = 1;
for i_fish = range_fish
    %% load fish
    ClusterIDs = [6,1];
    [cIX_load,gIX_load] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    cIX_abs = absIX(cIX_load);
    M_xyz_norm = CellXYZ_norm(cIX_abs,:);
    
    %%
    N = zeros(1,length(unique(gIX_load)));
    for i = 1:length(unique(gIX_load))
        IX = find(gIX_load==i);
        N(i) = length(IX);
    end
    %% exclude clusters with cell#>500
    temp = unique(gIX_load);
    U = temp(N<500);
    numClus = length(U);
    D = zeros(numClus,numClus);
    % set diag = realmax
    
    %%
    tic
    for i_clus_ref = 1:numClus
%         if mod(i_clus_ref,10)==0
%             disp(['i_clus_ref = ' num2str(i_clus_ref)]);
%         end
        
        clusID_ref = U(i_clus_ref);
        IX = find(gIX_load == clusID_ref);
        XYZ_ref = M_xyz_norm(IX,:);
        %
        %     for i_testnum = 1:numFish,%(i_refnum+1):numFish,% only do ordered pairwise test (other half redundant)
        %         i_fish_test = range_fish(i_testnum);
        %
        %         if i_fish_test == i_fish_ref,
        %             continue;
        %         end

        % cycle through clusters in test fish
        for i_clus_test = i_clus_ref:numClus%1:numClus % 
            if i_clus_test == i_clus_ref
                D(i_clus_ref,i_clus_test) = 0;
            else
                clusID_test = U(i_clus_test);
                IX = find(gIX_load == clusID_test);
                XYZ_test = M_xyz_norm(IX,:);
                
                % compute distance score
                score = ClusterDistanceD12(XYZ_ref,XYZ_test);
                D(i_clus_ref,i_clus_test) = score;
                D(i_clus_test,i_clus_ref) = score;
            end
        end              
    end
    toc
    %% color and plot anat
    D_small = zeros(1,numClus);
    for i_clus_ref = 1:numClus
        m = D(i_clus_ref,:);
        m_srt = sort(m);
        D_small(i_clus_ref) = mean(m_srt(2:6));
    end
%     figure;
%     bar(D_small)

    %%    
% clr1 = [1,0,0];
%     clr1_ = [0.5,0.4,0.4];
% %     clr1_ = [0.7,0.5,0.5];
%     numC = 64;
%     clrmap = Make1DColormap([clr1_;clr1],numC);
%     clrmap = rand(64,3);
    %%
    
    [cIX_select,gIX_select] = SelectClusterRange(cIX_load,gIX_load,U);
    gIX_select = SqueezeGroupIX(gIX_select);
    
    clrIX = MapXto1Dcolormap(D_small,[1,6],64);
    gIX = gIX_select;
    for i = 1:numClus
        gIX(gIX_select==i) = clrIX(i);
    end

    cIX = cIX_select;
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
    %     I.clrmap = clrmap;
    [h,im_full] = DrawCellsOnAnat(I);
    
    % save figure
    close(h);
    IM_full{1,i_fish} = im_full;
    
    %% select more isolated clusters, for a separate average plot
    [cIX,gIX] = SelectClusterRange(cIX,gIX,11:63);

    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
    %     I.clrmap = clrmap;
    [h,im_full] = DrawCellsOnAnat(I);
    
    % save figure
    close(h);
    IM_full{2,i_fish} = im_full;
end

M_reg_name = {'anatD12rank','anatD12rank_isolated'};
%% save as tiff stack
for i_reg = 1:n_reg
    range_im = range_fish;%1:18;
    tiffdir = fullfile(outputDir,[M_reg_name{i_reg} '_allfish.tiff']);
    % tiffdir = fullfile(outputDir,'White_1reg_allfish.tiff');
    
    % display each plane and save as tif
    h = figure;
    for i_plane = range_im
        im = IM_full{i_reg,i_plane};
        image(im);axis equal; axis off
        drawnow;
        % save tiff
        if (i_plane == 1)
            imwrite(im, tiffdir, 'compression','none','writemode','overwrite')
        else
            imwrite(im, tiffdir, 'compression','none','writemode','append')
        end
        %     pause(0.2)
    end
    close(h)
end

%% save avrage plot
for i_reg = 1:n_reg
    range_im = range_fish;%[1:3,5:18];
    cellarray = IM_full(i_reg,range_im);
    
    k_scale = 0.6; % this changes for every reg pair...
    k_contrast = 1;
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    imwrite(im_avr, fullfile(outputDir,[M_reg_name{i_reg} '_avr.tiff']), 'compression','none','writemode','overwrite');
end
