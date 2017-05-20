% fig3B top x% stimlock
clear all; close all; clc

isHotmapNotSingleColor = 1;
topPrct = 10; % 2
bottomPrct = 90; % 95

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range = 1:18;
M_reg_name = {'top10_allfish','bottom10_allfish'};
n_reg = length(M_reg_name);

IM_full = cell(n_reg,18);
for i_fish = range
    ClusterIDs = GetClusterIDs('all');
    [~,stimrange] = GetStimRange('2',i_fish);
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    
    [~,~,M_score] = GetTrialAvrLongTrace(hfig,M);
    
    [score_sort, IX_sort] = sort(M_score);
    
    %% set colormap
    if isHotmapNotSingleColor
        numC = 64;
        clrmap = flipud(hot(numC));
    else
        clr1 = [1,0,0];
        clr1_ = [0,0,0];
        numC = 64;
        clrmap = Make1DColormap([clr1_;clr1],numC);
    end
    
    %% top %
    thres_cut = prctile(M_score,topPrct);
    thres_min = prctile(M_score,0.1);
    ix_end = find(score_sort<thres_cut,1,'last');
    
    IX = IX_sort(1:ix_end);
    cIX = cIX_load(IX);
    Xrange = [thres_min,thres_cut];
    
    gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
    
    
    % anat plot
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);
    close(h);
    % save
    IM_full{1,i_fish} = im_full;
    
    %% bottom %
    thres_cut = prctile(M_score,bottomPrct);
    thres_max = prctile(M_score,99.9);
    ix_start = find(score_sort>thres_cut,1,'first');
    
    IX = IX_sort(ix_start:end);
    cIX = cIX_load(IX);
    Xrange = [thres_cut,thres_max];
    numC = 64;
    gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
    clrmap = flipud(hot(numC));
    
    % anat plot
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);
    close(h);
    % save
    IM_full{2,i_fish} = im_full;
    
end

%% save as tiff stack
outputDir = GetOutputDataDir;
for i_set = 1:n_reg
    range_im = range;%M_fishrange{i_set};
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
    IM = IM_full(i_set,range_im);
    
    SaveImToTiffStack(IM,tiffdir);    
end

%% avr
M_scale = {0.5,0.5};
% M_scale = {1,0.7};
for i_set = 1:2;
    range_im = range;%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    % adjust params for visualization
    k_scale = M_scale{i_set}; 
    k_contrast = 1;% 1.2;
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
%     text(20,1190,['n=',num2str(length(range))],'color','w','fontsize',14');
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end