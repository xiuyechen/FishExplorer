% new procedure: generate motor seed/map in trialRes space
% find top 100 cells correlating to motor, screen with Rh4+Rh5(+Rh6??) masks

% corr sweep (0.5-0.7) to get top %2 cells, then kmeans, save; rank by stim-lock, save
clear all; close all; clc

% manual choices!!
islrRes = 0;
istRes = 0; 
istAvrMotor = 1;

range_fish = 8:17;%GetFishRange;%[1:3,5:18];
M_stimrange = GetStimRange('O');%('2');


% set thres fol colormap value lower bound
if istAvrMotor
    reg_thres = 0.25;
else
    reg_thres = 0.25;
end


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
M_reg_name{1} = 'motormap_OMR_notRes_motorseed2';
M_ClusterIDs{1} = [11,2];
M_fishrange{1} = range_fish;
M_stimrange{1} = GetStimRange('O');
n_reg = 1;

i_set = 1;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
global VAR;

setappdata(hfig,'isMotorseed',1); % default is 1
if istRes
    setappdata(hfig,'isTrialRes',1); % default is 0
end
%% run fish
M_thres_reg = zeros(3,18);
M_numTopCorr = zeros(1,18);
M_motorseedRegs = cell(1,18);
M_compareMotorCellNumber = zeros(2,18);

IM_full = cell(n_reg,18);

[M_fishrange_im,fishrange_load] = CheckIfLoadFish(M_fishrange,M_ClusterIDs,M_stimrange);
for i_fish = fishrange_load
    % load fish
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,M_ClusterIDs{i_set},M_stimrange{i_set});

    %% regression (in tRes), thresholding by % of cells (instead of corr thres)
    
    %     [~,Reg_tRes] = GetTrialAvrLongTrace(hfig,Reg);
    %     [~,M_0_tRes] = GetTrialAvrLongTrace(hfig,M_0);
    
    Reg = FindClustermeans(gIX_seed,M); % this is already tRes
    if islrRes
        Reg = Reg-repmat(mean(Reg),2,1);
    end
    if istAvrMotor
        [Reg,~] = GetTrialAvrLongTrace(hfig,Reg);
    end
    
    %% compute correlation and choose top % cells
    Corr = corr(Reg',M_0');
    
    % keep best regression only
    [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr);
    
    % top 2 %
    nCells_total = size(M_0,1);
    prct_const = 2;
    [CIX,RegThres] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_rows);
    cIX1 = CIX{1};
    cIX2 = CIX{2};
    
    % get map color
    clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres,1],64);
    clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres,1],64);
%     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres1,1],64);
%     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres2,1],64);
    
    cIX = [cIX1;cIX2];
%     if isempty(cIX)
%         M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish); %#ok<SAGROW>
%         continue;
%     end
    
    clrIX = [clrIX1;clrIX2];
    gIX_offset = [ones(size(cIX1));2*ones(size(cIX2))];
    gIX = clrIX+(gIX_offset-1)*64;
%     gIX = [clrIX1;64+clrIX2];
    numK = length(unique(gIX));
    
    %% pool stats
    M_thres_reg(1,i_fish) = RegThres{1};
    M_thres_reg(2,i_fish) = RegThres{2};
    M_compareMotorCellNumber(1,i_fish) = length(cIX1);
    M_compareMotorCellNumber(2,i_fish) = length(cIX2);
    
    %% make double colormap
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
    
    %% make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full] = DrawCellsOnAnat(I);    
    
    %% save figure
    close(h);
    IM_full{i_set,i_fish} = im_full;    

end

%% draw color bars - to save???
figure
ax = axes('Position',[0.75,0.8,0.05,0.15],'Units','normalized');
DrawCustomColorbar(clrmap1,[reg_thres,1],2,ax);

ax = axes('Position',[0.9,0.8,0.05,0.15],'Units','normalized');
DrawCustomColorbar(clrmap2,[reg_thres,1],2,ax);

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

for i_set = 1:n_reg;
    range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    % adjust params for visualization
    k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.1;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end

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
