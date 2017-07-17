
%%
tagname = 'CV_clus_pass5';
outputDir = GetOutputDataDir;

range_fish = 1:18;%GetFishRange;%[1:3,5:18];%8:17;%
IM_full = cell(1,18);

%% init
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
for i_fish = range_fish
    
    ClusterIDs = [4,1];
    [cIX1,gIX1,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    ClusterIDs = [5,1];
    [cIX2,gIX2] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2));
    
    
    [score,im1,cIX_int,gIX_A,gIX_B] = HungarianCV(cIX1,cIX2,gIX1,gIX2,0);
    
    d = diag(im1);
    
    range = find(d>=5)';
    [cIX,gIX] = SelectClusterRange(cIX1,gIX1,range);
    
    gIX = SqueezeGroupIX(gIX);
    
    %% rank by stim-lock
    isStimLockRanking = 1;
    if isStimLockRanking
        M = UpdateIndices_Manual(hfig,cIX,gIX);
        C = FindClustermeans(gIX,M);
        [~,~,H] = GetTrialAvrLongTrace(hfig,C);
        [gIX,rankscore] = SortGroupIXbyScore(H,gIX);
    else
        M = UpdateIndices_Manual(hfig,cIX,gIX);
        C = FindClustermeans(gIX,M);
        numU = max(gIX);
        [gIX,rankscore] = RankByMotorReg_Direct(hfig,gIX,numU,C,1);
    end
    
    %% draw anat
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
    [h,im_full] = DrawCellsOnAnat(I);
    close(h);
    IM_full{i_fish} = im_full;
    
end

%% save as tiff stack
range_im = range_fish;%M_fishrange_im{i_set};
tiffdir = fullfile(outputDir,[tagname,'_allfish.tiff']);
IM = IM_full(range_im);

SaveImToTiffStack(IM,tiffdir);

%% Average Plot

range_im = range_fish;%M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
cellarray = IM_full(range_im);

% adjust params for visualization
k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
k_contrast = 1.1;%M_k_contrast{i_set};

[h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);

tiffdir = fullfile(outputDir,[tagname,'_avr.tiff']);
imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
