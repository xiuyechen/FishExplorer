% This is a batch script for the multi-stim t-test.

clear all; close all; clc

%% Setup
% saveFigFlag = 1;
outputDir = GetOutputDataDir;

M_pair_range = {[3,2],[9,8]}; % [3,2] for PT L vs R; [9,8] for OMR L vs R
M_reg_name = {'PT_ttest_pVal_-3','OMR_ttest_pVal_-3','PT_AND_OMR_ttest_pVal_-3','PT_not_OMR_ttest','OMR_not_PT_ttest'};
% M_pair_stimrange = {[1,2],[1,2]};
n_pair = 2;%length(M_pair_range);
 

range_fish = 8:18;%GetFishRange;%1:18;%
IM = cell(n_pair+3,max(range_fish));
IM_full = cell(n_pair+3,max(range_fish));

thres_prct = 2;
% thres_ttest = 0.001;

%% Init
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Main

M_hb_count = zeros(18,10);

for i_fish = range_fish       
    
    ClusterIDs = GetClusterIDs('all');
    stimrange = [1,2];
    [cIX_load,gIX_load,M,stim,behavior] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);

    % init
    M_cIX = cell(1,n_pair);
    
    %% loop through the pairs
    for i_set = 1:n_pair
        
        %%
        stim = getappdata(hfig,'stim');
        fishset = getappdata(hfig,'fishset');
        tteststimrange = M_pair_range{i_set};
        
        [stimStateBinary_all, ~] = GetStimStateBinary(stim,fishset);
        stimStateBinary = stimStateBinary_all(tteststimrange,:);
        if isempty(stimStateBinary),
            return;
        end
        
        %% Misha's method of averaging over 5 frames;
        if 1,
            samples = cell(1,length(tteststimrange));
            for i = 1:length(tteststimrange),
                IX = find(stimStateBinary(i,:));
                len = 5*floor(length(IX)/5);
                M_ = M(:,IX(1:len));
                temp = reshape(M_,size(M,1),5,[]);
                samples{i} = squeeze(mean(temp,2));
            end
            
            nCells = size(M,1);
            pvals = zeros(1,nCells);
            signs = zeros(1,nCells);
            for i_cell = 1:nCells,
                [~, p] = ttest2(samples{1}(i_cell,:),samples{2}(i_cell,:));
                pvals(i_cell) = p;
                signs(i_cell) = sign(mean(samples{1}(i_cell,:)) -  mean(samples{2}(i_cell,:)));
            end
        else
            %% Linear Discrimination Analysis with HotellingT2 test
            disp('t-test with HotellingT2...');
            
            thres_ttest = 0.001;
            % reshape
            sampling_interval = 5; % (down-sampling the time-dimension)
            M_3D = cell(1,length(tteststimrange));
            nCells = size(M,1);
            for i = 1:length(tteststimrange),
                diffStim = horzcat(stimStateBinary(i,1)==1,diff(stimStateBinary(i,:)));
                startIX = find(diffStim==1);
                stopIX = find(diffStim==-1);
                nBlocks = min(length(startIX),length(stopIX));
                startIX = startIX(1:nBlocks);
                stopIX = stopIX(1:nBlocks);
                
                intervals = stopIX-startIX;
                period = mode(intervals);
                IX = find(intervals==period);
                if length(IX)<length(startIX),
                    startIX = startIX(IX);
                    stopIX = stopIX(IX);
                end
                period_effective = length(1:sampling_interval:period);
                M_3D{i} = zeros(nCells,period_effective,nBlocks);
                for i_block = 1:nBlocks,
                    M_3D{i}(:,:,i_block) = M(:,startIX(i_block):sampling_interval:stopIX(i_block)-1);
                end
            end
            
            % Test with 'HotellingT2' function
            pvals = zeros(1,nCells);
            tic
            for i_cell = 1:nCells,
                M_1 = squeeze(M_3D{1}(i_cell,:,:))'; % rows ~ observations; col ~ dimensions
                M_2 = squeeze(M_3D{2}(i_cell,:,:))'; % rows ~ observations; col ~ dimensions
                
                %         Miu1 = mean(M_1,1)';
                %         Covar1 = cov(M_1);
                %         Miu2 = mean(M_2,1)';
                %         Covar2 = cov(M_2);
                %         Sigma_inv = inv(0.5*(Covar1+Covar2));
                %         w = Sigma_inv * (Miu1 - Miu2);
                
                % format input matrix for HotellingT2
                X1 = horzcat(ones(size(M_1,1),1),M_1);
                X2 = horzcat(2*ones(size(M_2,1),1),M_2);
                X = vertcat(X1,X2);
                % direct call to branch-function of 'HotellingT2' package
                pvals(i_cell) = T2Hot2ihe_Direct(X,thres_ttest);
            end
            toc
        end
        
        %% histogram plot
        mean(pvals)
        std(pvals)
        % set lower bound
        pvals_ = pvals;
        pvals_(find(pvals==0)) = 10^-100;
%         figure;hold on;
%         xv = -100:1:0;
%         counts = hist(log(pvals_),xv);
%         bar(-100:1:0,counts);
                
        [B,IX] = sort(pvals_);
             
        % select top cells
        if 0
            thres_ttest = B(nCells_target);
            nCells_target = round(thres_prct/100 * nCells);   
        else
            thres_ttest = 0.001;%1.0000e-10;
            nCells_target = length(find(B<thres_ttest));
        end
%         plot([log(thres_ttest),log(thres_ttest)],[0,max(counts)],'r--');
        
        %% rank cells, and threshold if desired
%         [B,I] = sort(pvals_);
        cIX_sorted = cIX_load(IX);              
        cIX = cIX_sorted(1:nCells_target);            
        gIX_raw = log(B(1:nCells_target));
        signs_ = signs(IX(1:nCells_target));
        gIX_signed = gIX_raw.*signs_;

%         if thres_ttest>0,
%             ix = find(B<thres_ttest,1,'last');
%             cIX = cIX(1:ix);
%         end
        
        % map p-values to colormap
%         gIX = ones(size(cIX));
        gIX = MapXto1Dcolormap(gIX_signed,[-50,50],64);
        % end
        
        clr1 = [1,0,0];
        clr0 = [0.3,0.3,0.3];
        clr1_ = [0,1,1];
%         clr1_ = [0.5,0.4,0.4];
        %     clr1_ = [0.7,0.5,0.5];
        numC = 64;
        clrmap = Make1DColormap([clr1_;clr0;clr1],numC);
        
        %% right plot
        I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
        I.clrmap = clrmap;
        [h,im_full,im] = DrawCellsOnAnat(I);
        
        %% save figure
        close(h);
        IM{i_set,i_fish} = im;
        IM_full{i_set,i_fish} = im_full;
        
        M_cIX{i_set} = cIX;
        
        %% count cells: hindbrain vs not
        MASKs = getappdata(hfig,'MASKs');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        
        %     Msk_IDs = [94,219,220];% midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
        cIX_mb = ScreenCellsWithMasks(94,cIX,gIX,MASKs,CellXYZ_norm,absIX);
        cIX_hb = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
        
        M_hb_count(i_fish,(i_set-1)*2+1) = length(cIX_mb);
        M_hb_count(i_fish,(i_set-1)*2+2) = length(cIX_hb);
    
    end
    
    %% 1: intersection (PT & OMR)
    i_set = 3;
    
    [cIX,ix] = intersect(M_cIX{1},M_cIX{2});
    %     gIX = M_gIX{1}(ix);
    gIX = 64*ones(length(cIX),1);
    %     [set2,ix2] = setdiff(CIX{2},CIX{1});
    %     cIX = [set1;set2];
    %     gIX = [GIX{1}(ix1);GIX{2}(ix2)];
    
    % make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full,im] = DrawCellsOnAnat(I);
    
    %% save figure
    close(h);
    IM{i_set,i_fish} = im;
    IM_full{i_set,i_fish} = im_full;
    
    %% count cells: hindbrain vs not
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    
    %     Msk_IDs = [94,219,220];% midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
    cIX_mb = ScreenCellsWithMasks(94,cIX,gIX,MASKs,CellXYZ_norm,absIX);
    cIX_hb = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
    
    M_hb_count(i_fish,n_pair*2+1) = length(cIX_mb);
    M_hb_count(i_fish,n_pair*2+2) = length(cIX_hb);
    
    %% 2: PT but not OMR
    i_set = 4;
    
    [cIX,ix] = setdiff(M_cIX{1},M_cIX{2});
    
    %     gIX = M_gIX{1}(ix);
    gIX = 64*ones(length(cIX),1);
    %     [set2,ix2] = setdiff(CIX{2},CIX{1});
    %     cIX = [set1;set2];
    %     gIX = [GIX{1}(ix1);GIX{2}(ix2)];
    
    % make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full,im] = DrawCellsOnAnat(I);
    
    %% save figure
    close(h);
    IM{i_set,i_fish} = im;
    IM_full{i_set,i_fish} = im_full;
    
    %% count cells: hindbrain vs not
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    
    %     Msk_IDs = [94,219,220];% midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
    cIX_mb = ScreenCellsWithMasks(94,cIX,gIX,MASKs,CellXYZ_norm,absIX);
    cIX_hb = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
    
    M_hb_count(i_fish,n_pair*2+3) = length(cIX_mb);
    M_hb_count(i_fish,n_pair*2+4) = length(cIX_hb);
    
    %% 3: OMR but not PT
    i_set = 5;
    
    [cIX,ix] = setdiff(M_cIX{2},M_cIX{1});
    %     gIX = M_gIX{1}(ix);
    gIX = 64*ones(length(cIX),1);
    %     [set2,ix2] = setdiff(CIX{2},CIX{1});
    %     cIX = [set1;set2];
    %     gIX = [GIX{1}(ix1);GIX{2}(ix2)];
    
    % make figure
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
    [h,im_full,im] = DrawCellsOnAnat(I);
    
    %% save figure
    close(h);
    IM{i_set,i_fish} = im;
    IM_full{i_set,i_fish} = im_full;
    
    %% count cells: hindbrain vs not
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    
    %     Msk_IDs = [94,219,220];% midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
    cIX_mb = ScreenCellsWithMasks(94,cIX,gIX,MASKs,CellXYZ_norm,absIX);
    cIX_hb = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
    
    M_hb_count(i_fish,n_pair*2+5) = length(cIX_mb);
    M_hb_count(i_fish,n_pair*2+6) = length(cIX_hb);
    
end

%% save as tiff stack
for i_set = 1:5
    range_im = range_fish;%1:18;
    tiffdir = fullfile(outputDir,[M_reg_name{i_set} '_allfish.tiff']);
    
    % display each plane and save as tif
    h = figure;
    for i_plane = range_im
        im = IM_full{i_set,i_plane};
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
%% [for later] plot from tiff stack
% isPlotfromtiffstack = 0;
% if isPlotfromtiffstack
%     IM_full = cell(1,18);
%     for i = 1:18
%         im = double(imread(tiffdir,i))./255;
%         IM_full{i} = im(317:1236,1:621,:);
%     end
% end

%% save avrage plot
for i_set = 3:5
    %%
%     i_set = 3;
    range_im = range_fish;%[1:3,5:18];
    cellarray = IM_full(i_set,range_im);
    
    k_scale = 0.3; % this changes for every reg pair...
    k_contrast = 2;
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    imwrite(im_avr, fullfile(outputDir,[M_reg_name{i_set} '_avr.tiff']),...
        'compression','none','writemode','overwrite');
end

%% count cell numbers

figure('Position',[100,400,150,160]);hold on;

inc1 = 0.2;
inc = 0.15;

x = 1;
Y = M_hb_count(fishrange,5)./(M_hb_count(fishrange,1)+M_hb_count(fishrange,3));
y1 = Y;
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
    0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

x = 2;
Y = M_hb_count(fishrange,6)./(M_hb_count(fishrange,2)+M_hb_count(fishrange,4));
y2 = Y;
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
    0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

xlim([0.5,2.5])
% ylim([0,2])
set(gca,'XTickLabels',{'midbrain','hindbrain'},'XTickLabelRotation',45);
ylabel('fraction of bimodal cells')

[h,p] = ttest2(y1,y2)

%% count cell numbers: setdiff, intersect
figure('Position',[100,400,250,200]);hold on;

inc1 = 0.2;
inc = 0.15;

fishrange = 8:18;
x = 1;
Y = (M_hb_count(fishrange,1)+M_hb_count(fishrange,3)-M_hb_count(fishrange,5));
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

x = 2;
Y = M_hb_count(fishrange,5);
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

x = 3;
Y = (M_hb_count(fishrange,2)+M_hb_count(fishrange,4)-M_hb_count(fishrange,6));
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

x = 4;
Y = M_hb_count(fishrange,6);
scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);

xlim([0.5,4.5]);
% ylim([0,2])
set(gca,'XTick',1:4,'XTickLabels',{'mB: single stimulus','mB: multi-stim',...
    'hB: single stimulus','hB: multi-stim'},'XTickLabelRotation',45);
ylabel('# of cells')


%% count cell numbers: all 10 categories including setdiff, intersect

figure('Position',[100,400,350,160]);hold on;

inc1 = 0.2;
inc = 0.15;

M_col = [1,7,3,9,5,2,8,4,10,6];
nCol = length(M_col);
xlabels = {'PT','PT-not-OMR','OMR','OMR-not-PT','PT&OMR','PT',...
    'PT-not-OMR','OMR','OMR-not-PT','PT&OMR'};

for x = 1:nCol
    Y = M_hb_count(8:end,M_col(x));
    scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],...
        'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
    plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
    err = std(Y)/sqrt(length(Y));
    plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
    plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
    plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);
end
xlim([0.5,nCol+0.5])
% ylim([0,2])
set(gca,'XTick',1:nCol,'XTickLabels',xlabels,'XTickLabelRotation',45);
ylabel('# of cells')

%%

save(fullfile(outputDir,'M_hb_count.mat'),'M_hb_count')
