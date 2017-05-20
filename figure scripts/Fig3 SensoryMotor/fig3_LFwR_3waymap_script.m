% using the setdiff idea as in fig3_motormap_setdiff_script, we are
% plotting the forward swimming network as a contrast to (either) left/right
% directional networks. This is now relying on the original fictive traces.
% (To use motorseeds, we'd have to manually separate the left, right and
% forward seeds, not easy in most fish.)

clear all; close all; clc

%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
M_reg_name{1} = 'L_Fw_R_map_notmotorseed_highercontrast';

range_fish = GetFishRange;%[1:3,5:18];%8:17;%
% stimrange = 2;

M_fishrange_im{1} = range_fish;
n_reg = 1;

i_set = 1;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',0);

%% run fish

% M_thres_reg = zeros(3,18);
% M_numTopCorr = zeros(1,18);
% M_motorseedRegs = cell(1,18);
% M_compareMotorCellNumber = zeros(2,18);

IM_full = cell(n_reg,18);

%%
for i_fish = range_fish
    
    %% load
    ClusterIDs = [1,1];
    [~,~,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    %% regression
    Reg = behavior([1:3],:); % left/Fw/right
    
    Corr = corr(Reg',M_0');
    
    % keep best regression only
    [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
    
    % top 2 %
    nCells_total = size(M_0,1);
    prct_const = 2;
    [CIX,~] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_rows);
    
    % get map color
    reg_thres_low = 0.2;
    reg_thres_high = 0.5;
    clrIX1 = MapXto1Dcolormap(corr_max(CIX{1}),[reg_thres_low,reg_thres_high],64);
    clrIX2 = MapXto1Dcolormap(corr_max(CIX{2}),[reg_thres_low,reg_thres_high],64);
    clrIX3 = MapXto1Dcolormap(corr_max(CIX{3}),[reg_thres_low,reg_thres_high],64);
    %     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres1,1],64);
    %     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres2,1],64);
    
    cIX = [CIX{1};CIX{2};CIX{3}];
    %     if isempty(cIX)
    %         M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish); %#ok<SAGROW>
    %         continue;
    %     end
    
    clrIX = [clrIX1;clrIX2;clrIX3];
    gIX_offset = [ones(size(CIX{1}));2*ones(size(CIX{2}));3*ones(size(CIX{3}))];
    gIX = clrIX+(gIX_offset-1)*64;
    %     gIX = [clrIX1;64+clrIX2];
    numK = length(unique(gIX));
    
    
    %% make triple colormap
    % R/G/B for 3 channels
    clr1 = [1,0,0];
    clr1_ = [0.5,0.4,0.4];
    %     clr1_ = [0.7,0.5,0.5];
    clr2 = [0,1,0];
    clr2_ = [0.4,0.5,0.4];
    
    clr3 = [0,0,1];
    clr3_ = [0.4,0.4,0.5];
    %     clr2_ = [0.5,0.7,0.7];
    numC = 64;
    clrmap1 = Make1DColormap([clr1_;clr1],numC);
    clrmap2 = Make1DColormap([clr2_;clr2],numC);
    clrmap3 = Make1DColormap([clr3_;clr3],numC);
    clrmap = [clrmap1;clrmap2;clrmap3];
       
    %% make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);
        
    %% save figure
    close(h);
    IM_full{i_set,i_fish} = im_full;
    
end

%% draw color bars - to save???
% figure
% ax = axes('Position',[0.75,0.8,0.05,0.15],'Units','normalized');
% DrawCustomColorbar(clrmap1,[reg_thres,1],2,ax);
% 
% ax = axes('Position',[0.9,0.8,0.05,0.15],'Units','normalized');
% DrawCustomColorbar(clrmap2,[reg_thres,1],2,ax);

%% save as tiff stack
for i_set = 1:n_reg
    range_im = M_fishrange_im{i_set};
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
    IM = IM_full(i_set,range_im);
    
    SaveImToTiffStack(IM,tiffdir);
end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

% optional override::::
% M_fishrange_im{1} = [1,3,5:17];
% M_fishrange_im{1} = 8:17;% for OMR
for i_set = 1:n_reg;
    range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    % adjust params for visualization
    k_scale = 0.8;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.3;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end

%%
if 0
    %% [for later] plot from tiff stack
    isPlotfromtiffstack = 0;
    
    if isPlotfromtiffstack
        IM_full = cell(n_reg,18);
        for i_set = 1:n_reg
            %% get tiff-stack path
            tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_allfish.tiff']);
            
            %% load
            
            for i_fish = 1:18
                im = double(imread(tiffdir,i_fish))./255;
                IM_full{i_set,i_fish} = im(317:1236,1:621,:);
            end
        end
    end
end

%%
% range_fish excludes Fish 4
% M_compareMotorCellNumber(:,4) = NaN;
% figure;bar(M_compareMotorCellNumber')
