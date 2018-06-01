clear all; close all; clc


%% folder setup
isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;


ClusterIDs = [2,1]; % init; can overrride
prct_const = 2; % init; can overrride



%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',1);

%%
IM_scatter = cell(2,18); % L and R

IM_maps = cell(5,18);
% IM_Xonly = cell(n_reg,18);
% IM_Yonly = cell(n_reg,18);
% IM_XY = cell(n_reg,18);
% IM_X = cell(n_reg,18);
% IM_Y = cell(n_reg,18);

caseID = 5;
switch caseID
    case 1
        load(fullfile(outputDir,'4D_SM_betas.mat'));
        M_reg_name = {'2x2motormaps'};
        range_fish = GetFishRange;% init; can overrride
    case 2
        load(fullfile(outputDir,'4D_SM_stimrangePT_betas.mat'));
        M_reg_name = {'2x2motormaps_PT'};
        range_fish = 6:18;
    case 3        
        load(fullfile(outputDir,'4D_SM_stimrangeOMR_betas.mat'));
        M_reg_name = {'2x2motormaps_OMR'};
        range_fish = 8:18;
    case 4
        load(fullfile(outputDir,'4D_SM_stimrangelooming_betas.mat'));
        M_reg_name = {'2x2motormaps_looming'};
        range_fish = [9:15,17:18];
    case 5
        load(fullfile(outputDir,'4D_SM_stimrangeDF_betas.mat'));
        M_reg_name = {'2x2motormaps_DF'};
        range_fish = [12:15,17:18];
end
%% run fish
for i_fish = range_fish
    
    cIX_all = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,[],0);
    
    %% main loop for left/right motor
    M_pass = cell(5,2);
    
    for i_lr = 1:2
        %% scatter plot with 4D components
        % loaded betas for this fish for the combo data (including both stim)
        betas = Betas{i_lr,i_fish};
        % set up plot dimensions
        X = betas(:,1);%b3;%b1;
        Y = betas(:,2);
        %                 Y = sqrt(b2.^2+b3.^2);%b2;
        Xname = 'motor res.';
        Yname = 'motor avr.';
        
        numcell = length(X);
        
        A = X;
        topN = round(prct_const/100*numcell); % top x% cutoff
        [~,IX] = sort(A,'descend');
        thresA = A(IX(topN));
        
        B = Y;
        topN = round(prct_const/100*numcell);
        [~,IX] = sort(B,'descend');
        thresB = B(IX(topN));
        
        IX_passX = find(A>=thresA);
        IX_passY = find(B>=thresB);
        IX_passXonly = setdiff(IX_passX,IX_passY);
        IX_passYonly = setdiff(IX_passY,IX_passX);%find(B>=thresB);%
        IX_passXY = intersect(IX_passX,IX_passY);
        IX_pass = union(IX_passX,IX_passY);
        IX_fail = intersect(find(A<thresA),find(B<thresB));%find(A<thresA);
        
        % get min/max
        x0 = min(X(IX_pass));
        x1 = max(X(IX_pass));
        y0 = min(Y(IX_pass));
        y1 = max(Y(IX_pass));
        
        gIX_in = (1:length(X))';
        
        M_pass{1,i_lr} = IX_passX;
        M_pass{2,i_lr} = IX_passY;
        M_pass{3,i_lr} = IX_passXonly;
        M_pass{4,i_lr} = IX_passYonly;
        M_pass{5,i_lr} = IX_passXY;
        %         PassX_2{i_lr} = IX_passX;
        %         PassY_2{i_lr} = IX_passY;
        %         PassXonly_2{i_lr} = IX_passXonly;
        %         PassYonly_2{i_lr} = IX_passYonly;
        %         PassXY_2{i_lr} = IX_passXY;
        
        %% scatter plot
        clrX = [0.3,0.8,0];
        clrY = [0.9,0.2,0.9];
        clrXY = [0,0.3,0.9];
        clr_fail = [0.5,0.5,0.5];
        
        h = figure('Position',[500,100,300,250]); hold on
        scatter(X,Y,1,clr_fail,'filled')

        scatter(X(IX_passXonly),Y(IX_passXonly),1,clrX);%,'filled');
        scatter(X(IX_passYonly),Y(IX_passYonly),1,clrY);%,'filled');
        scatter(X(IX_passXY),Y(IX_passXY),1,clrXY);%[1,0.5,0.5]);%,'filled');
        
        plot([x0,x1],[thresB,thresB],'k--');
        plot([thresA,thresA],[y0,y1],'k--');
        
        xlabel(Xname);ylabel(Yname);
        axis equal
%         set(gca,'XTick',0:0.5:1,'YTick',-0.5:0.5:1);
        axis tight
        %% save bubble plot
        IM_scatter{i_lr,i_fish} = print('-RGBImage');
        close(h)
    end
    
    %% component anat map
    setappdata(hfig,'clrmap_name','hsv_old');
    for i_map = 1:5
        cIX12 = M_pass(i_map,:);
        cIX = vertcat(cIX12{1},cIX12{2});
        gIX = vertcat(ones(size(cIX12{1})),2*ones(size(cIX12{2})));
        I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
        [h,~,im] = DrawCellsOnAnat(I);
        close(h);
        IM_maps{i_map,i_fish} = im;
    end
    
end


%% save as tiff stack
M_lr = {'-L','-R'};

M_comp_names = {'passX','passY','passXonly','passYonly','passXY'};
n_reg = 1;
for i_set = 1:n_reg
    range_im = range_fish;%M_fishrange_im{i_set};
    
    for i_lr = 1:2
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},M_lr{i_lr},'_scatter_allfish.tiff']);
        IM = IM_scatter(i_lr,range_im);
        SaveImToTiffStack(IM,tiffdir);
    end
    
    for i_map = 1:5
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_',M_comp_names{i_map},'_anat_allfish.tiff']);
        IM = IM_maps(i_map,range_im);
        SaveImToTiffStack(IM,tiffdir);
    end
end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

for i_set = 1:n_reg
    range_im = range_fish;%M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    
    for i_map = 1:5
        %%
        IM = IM_maps(i_map,range_im);%IM_2(i_lr,range_im);
        
        % adjust params for visualization
        k_scale = 0.5;%1/1.5;%M_k_scale{i_set};
        k_contrast = 1;%M_k_contrast{i_set};
        
        [h_anat,im_avr] = AverageAnatPlot(IM,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_',M_comp_names{i_map},'_anat_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
end
