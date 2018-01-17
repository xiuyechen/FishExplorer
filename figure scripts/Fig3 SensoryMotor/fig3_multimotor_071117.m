clear all; close all; clc

%% folder setup
outputDir = GetOutputDataDir;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];

%% set params!
% isCellbased = true;
ClusterIDs = [2,1];
% ClusterIDs = [7,1];

%%
tscriptstart = tic;
n_reg = 3;
IM_full = cell(n_reg,18);

for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    %% load data for chosen stim range
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);

    %% Method: stimAvr + motor regs
%     if isCellbased
        gIX = (1:length(cIX_load))';
        Data = M;
%     else % cluster based
%         gIX = gIX_load;
%         C = FindClustermeans(gIX,M);
%         Data = C;
%     end

    i_lr = 2;
    
    Data_tAvr = GetTrialAvrLongTrace(hfig,Data);
    [motor_tAvr,motor_tRes] = GetTrialAvrLongTrace(hfig,behavior);
    stimcorr = DiagCorr(Data_tAvr',Data');
    motorcorr = corr(motor_tRes(i_lr,:)',Data');

    xlims = diag(corr(motor_tAvr',behavior'));
    ylims = diag(corr(motor_tRes',behavior'));
    
    % make figures
    cIX_in = cIX_load;
    gIX_in = gIX;
    xbound = xlims(i_lr);
    ybound = ylims(i_lr);
    
    %%
    %% make custom 2-D colormap
    grid = Make4color2Dcolormap;
    % grid = MakeDiagonal2Dcolormap;
    
    %% get new gIX with matching custom colormap 'cmap_U'
    if false % temp, trying this for cell-based betas
        thres_stim = prctile(stimcorr,90);
        thres_motor = prctile(motorcorr,90);
        IX_stim = find(stimcorr>thres_stim);
        IX_motor = find(motorcorr>thres_motor);
        stimcorr(IX_stim) = thres_stim;
        motorcorr(IX_motor) = thres_motor;
    end
        
    %% threshold
    topN = 1000;
    [~,IX] = sort(motorcorr,'descend');
    y0 = motorcorr(IX(topN));
    x0 = min(stimcorr(IX(1:topN)));
    x1 = max(stimcorr(IX(1:topN)));
    y1 = motorcorr(IX(1));
    
    IX_pass = find(motorcorr>=y0);
    IX_fail = find(motorcorr<y0);
    %     IX_pass = find(motorcorr./(stimcorr-x0) > tang);
    
    disp(['# of units pass thres: ',num2str(length(IX_pass))]);
    [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_pass);
    
%     topN = 100;
%     [~,IX] = sort(motorcorr,'descend');
%     y1 = mean(motorcorr(IX(1:topN)));
%     x1 = mean(stimcorr(IX(1:topN)));
%     [~,IX] = sort(stimcorr);
%     x0 = mean(stimcorr(IX(1:topN)));
%     tang = y1/(x1*1.1-x0*1.1);

    %% map data to colormap, and cluster sizes
    clrmap0 = MapXYto2Dcolormap(gIX_in,stimcorr,motorcorr,[x0,x1],[y0,y1],grid);
    clrmap = clrmap0;
    clrmap(IX_fail,:) = ones(length(IX_fail),3)*0.5;
%     clrX_max = max(0.4,xbound);
%     clrY_max = max(0.7,ybound);
%     clrmap = MapXYto2Dcolormap(gIX_in,stimcorr,motorcorr,[0,clrX_max],[0.3,clrY_max],grid);
    
    if length(clrmap)>1000
        U_size = ones(size(stimcorr));
    else % for actual clusters
        
        % get cluster size (number of cells in each cluster)
        gIX2 = SqueezeGroupIX(gIX_in);
        U = unique(gIX2);
        U_size = zeros(size(stimcorr));
        for i = 1:length(U)
            ix = find(gIX2 == U(i));
            U_size(i) = length(ix);
        end
    end

    %%
    h = figure('Position',[500,500,300,400]); hold on
    scatter(stimcorr(IX_pass),motorcorr(IX_pass),3,clrmap(IX_pass,:),'filled')
    xlabel('tAvr corr.');ylabel('tRes corr.');
        
    scatter(stimcorr(IX_fail),motorcorr(IX_fail),1,clrmap(IX_fail,:));
    plot([x0*0.8,x1*1.2],[y0,y0],'k--');
    
    xlim([0,1])
    ylim([-0.5,1])
    axis equal
    set(gca,'XTick',0:0.5:1,'YTick',-0.5:0.5:1);
    %% bubble plot
    h = figure('Position',[500,500,300,250]); hold on
    scatter(stimcorr,motorcorr,U_size,clrmap,'filled')    
        xlabel('tAvr corr.');ylabel('tRes corr.');
    
    if 1
        axis equal
        xlim([0,1]);
        ylim([-0.3,1]);
        %     set(gca,'YTick',-0.2:0.2:0.6);
    end    
    
    plot([x0,x1],[y0,y0],'k--');
%     scatter(x1,y1,5,'k','filled');
%     plot([x0*1.1,x1*1.1],[0,y1],'k--');
    
    IM_full{1,i_fish} = print('-RGBImage');
    close(h)
    % plot 'limits' from behavior trace
%     plot([xbound,xbound],[-0.3,1],'m:');
%     plot([-0.3,1],[ybound,ybound],'g:');
%     plot([-0.3,1],[thres_plot,thres_plot],'r');
        
    %% anat, all cells included in A0.5       
    I = LoadCurrentFishForAnatPlot(hfig,cIX_in,gIX_in,clrmap,[]);
    [h,im] = DrawCellsOnAnat(I);
    close(h)
    IM_full{2,i_fish} = im;   
    
    %% anat, pass thres    
    I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
    [h,im] = DrawCellsOnAnat(I);
    close(h)
    IM_full{3,i_fish} = im;           
    
end
toc(tscriptstart)

M_reg_name = {'bubbleplot_R_all_passy','multimotor_R_anat_all','multimotor_R_anat_all_passy'};
% M_reg_name = {'bubbleplot_L_all_passy','multimotor_L_anat_all','multimotor_L_anat_all_passy'};
% M_reg_name = {'bubbleplot_R_A0.5_passy','multimotor_R_anat_A0.5','multimotor_R_anat_A0.5_passy'};
% M_reg_name = {'bubbleplot_R_all','multimotor_R_anat_all','multimotor_R_anat_all_passangle'};
% M_reg_name = {'bubbleplot_R','multimotor_R_anat_allA0.5','multimotor_R_anat_passangle'};
%% save as tiff stack
for i_set = 1:n_reg
    range_im = range_fish;%M_fishrange_im{i_set};
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
    range_im = [1:3,5:7];%range_fish;% M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    %% adjust params for visualization
    k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.1;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    %%
    tiffdir = fullfile(outputDir,[M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end

