% fig3B top 1% stimlock
%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range = 1:18;%GetFishRange;
IM_full = cell(2,18);
for i_fish = 9:18%range
    ClusterIDs = GetClusterIDs('all');
    %     stimrange = 1;
    % [cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    [~,~,M_score] = GetTrialAvrLongTrace(hfig,M);
    
    [score_sort, IX_sort] = sort(M_score);
    
    %% top 2%
    thres_cut = prctile(M_score,2);
    thres_min = prctile(M_score,0.1);
    ix_end = find(score_sort<thres_cut,1,'last');
    
    IX = IX_sort(1:ix_end);
    cIX = cIX_load(IX);
    Xrange = [thres_min,thres_cut];
    numC = 64;
    gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
    clrmap = flipud(hot(numC));
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);
    close(h);
    IM_full{1,i_fish} = im_full;
    
    %% bottom 5%
    thres_cut = prctile(M_score,95);
    thres_max = prctile(M_score,99.9);
    ix_start = find(score_sort>thres_cut,1,'first');
    
    IX = IX_sort(ix_start:end);
    cIX = cIX_load(IX);
    Xrange = [thres_cut,thres_max];
    numC = 64;
    gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
    clrmap = flipud(hot(numC));
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);
    close(h);
    IM_full{2,i_fish} = im_full;
end

%% save as tiff stack
for i_set = 1:n_reg
    range_im = M_fishrange{i_set};
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
    IM = IM_full(i_set,range_im);
    
    SaveImToTiffStack(IM,tiffdir);    
end

%%
M_reg_name = {'top2_allfish','bottom5_allfish'};
for i_set = 1:2;
    range_im = range;%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    % adjust params for visualization
    k_scale = 1; 
    k_contrast = 1;% 1.2;
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end