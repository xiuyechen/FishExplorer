% this code follows MultiSensoryConvergence\step1_savemats\Batch_PTvsOMR_regbased_save_intersection.m
% plots the saved intersections on anat map, color-coded for same or
% opposite direction (red,cyan=congruent; yellow,purple=incongruent).

% e.g. load PTintDF_regbased_11fish_sweepthres_cong.mat

clearvars; close all; clc

%% init
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

% setappdata(hfig,'clrmap_name','hsv_old');

%%
for caseflag = 1:6
    switch caseflag
        
        case 1 % PT & OMR
            M_reg_name{1} = 'PTintOMR_regbased';
            M_reg_range = {[3,2],[9,8]};
            M_stimrange = {1,2};
            range_fish = [8,9,11:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
            
        case 2 % PT & Looming
            M_reg_name{1} = 'PTintLm_regbased';
            M_reg_range = {[3,2],[11,12]};
            M_stimrange = {1,5};
            range_fish = [9:15,17:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
            
        case 3 % PT & DF
            M_reg_name{1} = 'PTintDF_regbased';
            M_reg_range = {[3,2],[1,4]}; % [3,1],[2,1] "congruent", [3,4],[2,4] "incongruent" % 3=PT-L, 1=DF
            M_stimrange = {1,3};
            range_fish = [1:5,12:15,17:18];%[12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]}; % PT&black, i.e. R-off, L-off
            M_incong_pairs = {[1,2],[2,2]}; % PT&white i.e. L-on, R-on
            
        case 4 % OMR & Looming
            M_reg_name{1} = 'OMRintLm_regbased';
            M_reg_range = {[9,8],[11,12]};
            M_stimrange = {2,5};
            range_fish = [9,11:15,17:18];
            
            M_cong_pairs = {[1,1],[2,2]};
            M_incong_pairs = {[1,2],[2,1]};
            
        case 5 % OMR & DF
            M_reg_name{1} = 'OMRintDF_regbased';
            M_reg_range = {[9,8],[11,12]};
            M_stimrange = {2,3};
            range_fish = [12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]};
            M_incong_pairs = {[1,2],[2,2]};
            
        case 6 % Looming & DF
            M_reg_name{1} = 'LmintDF_regbased';
            M_reg_range = {[3,2],[1,4]};
            M_stimrange = {5,3};
            range_fish = [12:15,17:18];
            
            M_cong_pairs = {[1,1],[2,1]};
            M_incong_pairs = {[1,2],[2,2]};
            
    end
    
    outputDir = GetOutputDataDir;
    load(fullfile(outputDir,[M_reg_name{1},'_sweepthres_cong.mat']),'Intersect_cIX','Intersect_gIX');
    
    %%
    i_prct_count = 5;
    reg_name{1} = [M_reg_name{1},'_',num2str(i_prct_count),'%'];
    
    n_reg = 1;
    i_set = 1;
    
    IM_full = cell(1,18);
    
    % make custom colormap
    clrmap0 = hsv(4);
    clrmap = clrmap0([1,3,2,4],:);
    clrmap(3,:) = [1,1,0];
    
    for i_fish = range_fish
        LoadFullFish(hfig,i_fish,0);
        
        cIX_int = Intersect_cIX{i_fish,i_prct_count};
        gIX_int = Intersect_gIX{i_fish,i_prct_count};
        
        % draw hindbrain Rh1/2/3 mask
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        setappdata(hfig,'isShowMasks',1);
        setappdata(hfig,'Msk_IDs',[219,220,221]);
        setappdata(hfig,'isShowMskOutline',1);
        
        % make figure
        I = LoadCurrentFishForAnatPlot(hfig,cIX_int,gIX_int,clrmap);
        [h,~,im_full1] = DrawCellsOnAnat(I);
        close(h);
        IM_full{i_set,i_fish} = im_full1;
        
    end
    
    %% save as tiff stack
    for i_set = 1:n_reg
        range_im = range_fish;%M_fishrange_im{i_set};
        tiffdir = fullfile(outputDir,[reg_name{i_set},'_allfish.tiff']);
        IM = IM_full(i_set,range_im);
        
        SaveImToTiffStack(IM,tiffdir);
    end
    
    %% Average Plot
    for i_set = 1:n_reg
        range_im = range_fish;% M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
        cellarray = IM_full(i_set,range_im);
        
        %% adjust params for visualization
        %         k_scale = 0.5;%1;%1/1.5;%M_k_scale{i_set};
        %         k_contrast = 0.9;%M_k_contrast{i_set};
        %         [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray);
        close(h_anat);
        
        %%
        tiffdir = fullfile(outputDir,[reg_name{i_set},'_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
end

%% Count intersection cell number, plot ratio into triangular matrix map
Ratio_intersect = zeros(4,4,18);
Ratio_cong = zeros(4,4,18);
Ratio_incong = zeros(4,4,18);

outputDir = GetOutputDataDir;

% i_prct_count = 5;
for i_fish = 1:18
    LoadFullFish(hfig,i_fish,0);
    absIX = getappdata(hfig,'absIX');
    nCells_total = length(absIX);
    nCells_target = round(i_prct_count/100 * nCells_total);
    
    %%
    for caseflag = 1:6
        switch caseflag
            case 1 % PT & OMR
                M_reg_name{1} = 'PTintOMR_regbased';
                ii = 1;
                jj = 2;
            case 2 % PT & Looming
                M_reg_name{1} = 'PTintLm_regbased';
                ii = 1;
                jj = 3;
            case 3 % PT & DF
                M_reg_name{1} = 'PTintDF_regbased';
                ii = 1;
                jj = 4;
            case 4 % OMR & Looming
                M_reg_name{1} = 'OMRintLm_regbased';
                ii = 2;
                jj = 3;
            case 5 % OMR & DF
                M_reg_name{1} = 'OMRintDF_regbased';
                ii = 2;
                jj = 4;
            case 6 % Looming & DF
                M_reg_name{1} = 'LmintDF_regbased';
                ii = 3;
                jj = 4;
        end
        load(fullfile(outputDir,[M_reg_name{1},'_sweepthres_cong.mat']),'Intersect_cIX','Intersect_gIX');
        cIX_int = Intersect_cIX{i_fish,i_prct_count};
        gIX_int = Intersect_gIX{i_fish,i_prct_count};
        
        if ~isempty(cIX_int)
            Ratio_intersect(ii,jj,i_fish) = length(cIX_int)/nCells_target;
            IX = find(gIX_int<=2); % 1 or 2
            Ratio_cong(ii,jj,i_fish) = length(IX)/nCells_target;
            IX = find(gIX_int>=3); % 3 or 4
            Ratio_incong(ii,jj,i_fish) = length(IX)/nCells_target;
        end
    end
end
%% visuals of the above
labels1 = {'phT','OMR','loom'};
labels2 = {'OMR','loom','DF'};

M_data = {Ratio_intersect,Ratio_cong,Ratio_incong};
for i_plot = 1:3
    Ratio = M_data{i_plot};
    Scores = nanmean(Ratio,3);
    Scores(isnan(Scores))=0;
    Scores = Scores(1:3,2:4);
    
    figure('Position',[100,400,250,180]);
    CorrPlot(Scores,1,labels1);
    
    set(gca,'XAxisLocation','top');
%     set(gca,'YAxisLocation','right');
    set(gca,'XTickLabel',labels2,'XTickLabelRotation',45)
    set(gca,'YTickLabelRotation',45);
%     box off
end

%% quantify by ZBrain masks
Ratio_intersect = zeros(4,4,18);
Ratio_cong = zeros(4,4,18);
Ratio_incong = zeros(4,4,18);

Msk_IDs = [275,1,94,219:225]; % forebrain 275; dien 1; midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
nMasks = length(Msk_IDs);
M_anat_count = zeros(18,nMasks,6);

% i_prct_count = 5;
for i_fish = 1:18
    LoadFullFish(hfig,i_fish,0);
    absIX = getappdata(hfig,'absIX');
    nCells_total = length(absIX);
    nCells_target = round(i_prct_count/100 * nCells_total);
    
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    
    %%
    for caseflag = 1:6
        switch caseflag
            case 1 % PT & OMR
                M_reg_name{1} = 'PTintOMR_regbased';
                ii = 1;
                jj = 2;
            case 2 % PT & Looming
                M_reg_name{1} = 'PTintLm_regbased';
                ii = 1;
                jj = 3;
            case 3 % PT & DF
                M_reg_name{1} = 'PTintDF_regbased';
                ii = 1;
                jj = 4;
            case 4 % OMR & Looming
                M_reg_name{1} = 'OMRintLm_regbased';
                ii = 2;
                jj = 3;
            case 5 % OMR & DF
                M_reg_name{1} = 'OMRintDF_regbased';
                ii = 2;
                jj = 4;
            case 6 % Looming & DF
                M_reg_name{1} = 'LmintDF_regbased';
                ii = 3;
                jj = 4;
        end
        load(fullfile(outputDir,[M_reg_name{1},'_sweepthres_cong.mat']),'Intersect_cIX','Intersect_gIX');
        cIX_int = Intersect_cIX{i_fish,i_prct_count};
        gIX_int = Intersect_gIX{i_fish,i_prct_count};
        
        if ~isempty(cIX_int)
            % count by anat masks
            for i_mask = 1:nMasks
                cIX = ScreenCellsWithMasks(Msk_IDs(i_mask),cIX_int,gIX_int,MASKs,CellXYZ_norm,absIX);
                M_anat_count(i_fish,i_mask,caseflag) = length(cIX)/length(cIX_int);
            end
            
            % intersection cell number ratios
            Ratio_intersect(ii,jj,i_fish) = length(cIX_int)/nCells_target;
            IX = find(gIX_int<=2); % 1 or 2
            Ratio_cong(ii,jj,i_fish) = length(IX)/nCells_target;
            IX = find(gIX_int>=3); % 3 or 4
            Ratio_incong(ii,jj,i_fish) = length(IX)/nCells_target;
        end
    end
end

%% multi-fish bar plot of the above
figure('Position',[100,400,450,860]);hold on;

inc1 = 0.2;
inc = 0.15;

M_cond = {'phT & OMR','phT & loom','phT & DF','OMR & loom','OMR & DF','loom & DF'};
for i_cond = 1:6
    subplot(6,1,i_cond);hold on;
    title(M_cond{i_cond});
    
    % connect the x columns with grey lines for each fish
    for i_fishcount = 1:18
        Y = squeeze(M_anat_count(i_fishcount,:,i_cond))*100;
        plot(1:nMasks,Y,'color',[0.7,0.7,0.7],'linewidth',0.5)
    end
    
    % plot the data points as dots, with avr and SEM
    for x = 1:nMasks
        %%
        Y = M_anat_count(:,x,i_cond)*100;
        Y = Y(find(Y));
        
        %     scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
        scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
            0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
        err = std(Y)/sqrt(length(Y));
        plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
        plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
        plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);
    end
    
    xlim([0.5,nMasks+0.5])
    ylim([0,75]);
    labels = {'Telen.','Dien.','midB','Rh1','Rh2',...
        'Rh3','Rh4','Rh5','Rh6','Rh7',};
    % labels = {'midbrain','hindbrain Rh1,2','hindbrain Rh3+'};
    set(gca,'XTick',1:nMasks,'XTickLabels',labels,'XTickLabelRotation',45);
    ylabel('% of int cells')
    
end