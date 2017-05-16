%% Setup
% folder setup
% saveFigFlag = 1;

outputDir = GetOutputDataDir;
% saveDir = [];
% saveDir{1} = fullfile(outputDir,'fig2_0422','PT_LR_stimrangePT');
% setDir(saveDir{1}); % make folder if doesn't exist
% saveDir{2} = fullfile(outputDir,'fig2_0422','LonRon_LR_stimrangePT');
% setDir(saveDir{2}); % make folder if doesn't exist
% saveDir{3} = fullfile(outputDir,'fig2_0422','motor_LR_stimrangePT');
% setDir(saveDir{3}); % make folder if doesn't exist

% params

reg_range = 4;
reg_thres = 0.5;
n_reg = 1;
% M_stimmotorflag = [1,1,0]; % 1 for stim and 0 for motor
% M_reg_range = {[3,2],[6,5],[4,4],[1,3]}; % [2,3] for PT, [8,9] for OMR for fish8; left/right pairs
% M_reg_thres = {0.5,0.5,0.5};
% n_reg = length(M_reg_thres);

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range = 1:18;%GetFishRange;
IM = cell(n_reg,max(range));
IM_full = cell(n_reg,max(range));
%%
i_set = 1;
for i_fish = range
    ClusterIDs = GetClusterIDs('all');
    stimrange = 1;
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    
    %% regression
    setappdata(hfig,'isStimAvr',0);
    [M_0,M,behavior,stim] = UpdateTimeIndex(hfig);
    
    fishset = getappdata(hfig,'fishset');
    [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
    
    Reg = regressors(reg_range,:);
    R = corr(Reg',M_0');
    
    
    cIX = (find(R>reg_thres))';
    wIX = R(cIX); % weight, i.e. corr.coeff
    
    
    %% set up custom colormap: hot (heatmap)
%     clrbins = 0.5:0.05:1;
%     clrmap0 = hot(round(1.5*length(clrbins)));
%     clrmap = clrmap0(end-length(clrbins):end,:);
%     
%     gIX = interp1(clrbins,1:length(clrbins),wIX,'nearest');
%         
    %% single color gradient (red)
    clr1 = [1,0,0];
    clr1_ = [0.7,0.5,0.5];
    numC = 64;
    clrmap = Make1DColormap([clr1_;clr1],numC);
        
    gIX = MapXto1Dcolormap(wIX,[reg_thres,1],64);
    
    %% right plot
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
    I.clrmap = clrmap;
    [h,im_full,im] = DrawCellsOnAnat(I);
    
    %% save figure
    close(h);
    IM{i_set,i_fish} = im;
    IM_full{i_set,i_fish} = im_full;
    %         savefolder = fullfile(saveDir{i_set},['regthres' num2str(reg_thres)]);
    %         setDir(savefolder);
    %         figName = ['Fish' num2str(i_fish)];
    % %         SaveFigureHelper(saveFigFlag, savefolder, figName);
    
end

%% save as tiff stack
i_reg = 1;
range_im = 1:18;
tiffdir = fullfile(outputDir,'White_1reg_allfish.tiff');

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

%% [for later] plot from tiff stack
isPlotfromtiffstack = 0;
if isPlotfromtiffstack
    IM_full = cell(1,18);
    for i = 1:18
        im = double(imread(tiffdir,i))./255;
        IM_full{i} = im(317:1236,1:621,:);
    end
end

%%
i_reg = 1;
range_im = 1:18;%,5:7];%[1:3,5:18];
cellarray = IM_full(i_reg,range_im);

k_scale = 2; % this changes for every reg pair...
k_contrast = 1.2;

[h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
imwrite(im_avr, fullfile(outputDir,'White_1reg_allfish_avr.tiff'), 'compression','none','writemode','overwrite');
%%
% anat_YX = getappdata(hfig,'anat_yx_norm');
% z_range_ventral = 40;
% z_range = z_range_ventral:345; % max 138*k_zres
% y_range = 81:1000;%81:990;% 81:1104; % max 1406
% dimv_yx3 = size(anat_YX(y_range,:,:));
% dimv_zx3 = [length(z_range),size(anat_YX,2),3];
% im_full = IM_full{1};
% im = IM_full{1}(dimv_zx3(1)+11:end,1:dimv_yx3(2),:);

%% function IM_crop = CropFullAnat(IM_full)
IM_crop = cell(size(IM_full));
for i = 1:length(IM_full)
    IM_crop{i} = IM_full{1}(317:1221,1:621,:);%(317:1236,1:621,:);
end

cellarray = IM_crop(1,1:4);
[h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);


