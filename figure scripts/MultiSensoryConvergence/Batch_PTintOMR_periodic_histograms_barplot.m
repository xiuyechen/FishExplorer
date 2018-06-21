%%
clear all; close all; clc

%% load

caseflag = 2;
switch caseflag 
    case 1
        outputDir = GetOutputDataDir;
        load(fullfile(outputDir,'PTintOMR_regbased_sweepthres'),'Intersect_cIX');
        Intersect_fixedprct = Intersect_cIX(:,5);
        % load(fullfile(outputDir,'PTintOMR_sweepthres'),'Intersect_cIX');
        % PTintOMR = Intersect_cIX(:,7);
        
        load(fullfile(outputDir,'4D_SM_stimrangePTOMR_minmax_betas.mat'));
        % load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));
        
        range_fish = 8:18;
        i_fish_example = 8;

    case 2
        outputDir = GetOutputDataDir;
        load(fullfile(outputDir,'OMRintLm_regbased_sweepthres'),'Intersect_cIX');
        Intersect_fixedprct = Intersect_cIX(:,5);
        % load(fullfile(outputDir,'PTintOMR_sweepthres'),'Intersect_cIX');
        % PTintOMR = Intersect_cIX(:,7);
        
        load(fullfile(outputDir,'4D_SM_stimrangeOMRlooming_minmax_betas.mat'));
        % load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));
        
        range_fish = [9:15,17:18];
        i_fish_example = 9;
        
end

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
XY_LR = cell(2,2,18);
XY_LR_IXpass = cell(2,2,18);
for i_fish = range_fish
    cIX_int = Intersect_fixedprct{i_fish};
    %%
    for i_ax = 1:2
        for i_lr = 1:2
            % get betas for this fish
            betas = Betas{i_lr,i_fish};
            b1 = betas(:,1);
            b2 = betas(:,2);
            b3 = betas(:,3);
            
            %     X = b1;%b3;%b1;
            %     Y = sqrt(b2.^2+b3.^2);%b2;
            %     Xname = 'motor only (b1)';
            %     Yname = 'periodic';
            
            numcell = length(b1);
            
            if i_ax==1
                A = b1;
            else
                A = betas(:,5);%b2;;
%                 A = sqrt(b2.^2+b3.^2);%b2;;
            end
            
            topN = length(cIX_int);%round(0.01*numcell); % top 5% cutoff
            [~,IX] = sort(A,'descend');
            thresA = A(IX(topN));
            IX_pass = find(A>=thresA);
            
            XY_LR{i_ax,i_lr,i_fish} = A;
            XY_LR_IXpass{i_ax,i_lr,i_fish} = IX_pass;
        end
    end
end

%% histogram: Distribution of convergence cells compared to all cells

% clr_XY = [0.5,1,0.5; 0.5,0.5,1];
clr_XY = [0.3,0.8,0.2; 0.1,0.3,1];


xbins = 0:0.05:1;
numbins = length(xbins)-1;
offset = (xbins(2)-xbins(1))/2;

% pool for histogram
N_XY_pass = zeros(length(range_fish)*2,numbins,2);
N_XY_int = zeros(length(range_fish)*2,numbins,2);
N_XY_full = zeros(length(range_fish)*2,numbins,2);

% pool for bar plot
s_XY_pass = zeros(length(range_fish)*2,2,2); % mean&sem, i_ax
s_XY_int = zeros(length(range_fish)*2,2,2);
s_XY_full = zeros(length(range_fish)*2,2,2);

for i_fishcount = 1;%1:length(range_fish)
    i_fish = range_fish(i_fishcount);
    cIX_int = Intersect_fixedprct{i_fish};
    %%
    for i_ax = 1:2
        n_pass = cell(1,2);
        n_int = cell(1,2);
        n_full = cell(1,2);
        for i_lr = 1:2
            %% unload data
            IX_pass = XY_LR_IXpass{i_ax,i_lr,i_fish};
            y_full = XY_LR{i_ax,i_lr,i_fish};
            y_pass = y_full(IX_pass);
            y_int = y_full(cIX_int);
            
            %% normalize with min/max
            y0 = min(y_full);
            y1 = max(y_full);
            y_pass_n = (y_pass-y0)/(y1-y0);
            y_int_n = (y_int-y0)/(y1-y0);
            y_full_n = (y_full-y0)/(y1-y0);
            
            %% plot overlapping histogram
            if i_fish==8 && i_lr==2
                % plot single examples manually
                
                if i_ax==1
                    figure('Position',[500,500,250,100]);hold on;
                    h1 = histogram(y_pass_n);
                    h2 = histogram(y_int_n);
                    xlim([0,1])
                    xlabel('% value')
                    ylabel('a.u.')
                    set(gca,'XTick',[0,1]);
                    set(gca,'YTick',[]);
                else
                    figure('Position',[500,500,100,250]);hold on;
                    h1 = histogram(y_pass_n,'Orientation','horizontal');
                    h2 = histogram(y_int_n,'Orientation','horizontal');    
                    ylim([0,1]);
                    xlim([0,0.3])
                    ylabel('% value')
                    xlabel('a.u.')
                    set(gca,'YTick',[0,1]);
                    set(gca,'XTick',[]);
                end
                
                h1.Normalization = 'probability';
                h1.BinWidth = 0.025;
                h2.Normalization = 'probability';
                h2.BinWidth = 0.025;
                
                h1.FaceColor = [1,1,1];%clr_XY(i_ax,:);
                h1.EdgeColor = clr_XY(i_ax,:);
                h2.FaceColor = [1,0.5,0.5];
                h2.EdgeColor = [1,0.5,0.5];
            end
            
            %% counts for multifish bar plot
            s_XY_pass((i_fishcount-1)*2+i_lr,1,i_ax) = mean(y_pass_n);
            s_XY_int((i_fishcount-1)*2+i_lr,1,i_ax) = mean(y_int_n);
            s_XY_full((i_fishcount-1)*2+i_lr,1,i_ax) = mean(y_full_n);
            
            s_XY_pass((i_fishcount-1)*2+i_lr,2,i_ax) = std(y_pass_n);
            s_XY_int((i_fishcount-1)*2+i_lr,2,i_ax) = std(y_int_n);
            s_XY_full((i_fishcount-1)*2+i_lr,1,i_ax) = std(y_full_n);
            
            %% counts for multifish histogram
            [n_pass{i_lr},edges] = histcounts(y_pass_n,xbins);
            [n_int{i_lr},edges] = histcounts(y_int_n,xbins);%+offset);
            [n_full{i_lr},edges] = histcounts(y_full_n,xbins);
            
            N_XY_pass((i_fishcount-1)*2+i_lr,:,i_ax) = n_pass{i_lr};
            N_XY_int((i_fishcount-1)*2+i_lr,:,i_ax) = n_int{i_lr}/2;
            N_XY_full((i_fishcount-1)*2+i_lr,:,i_ax) = n_full{i_lr}/2;
        end
    end
end

%% to make corresponding 2-D scatter plot related to the above histogram, use code in Batch_PTintOMR_periodic_scatterplot_anatmap.m

flag = true;
if flag
    %% example
    i_fish = i_fish_example; 
    i_lr=2;
    cIX_int = Intersect_fixedprct{i_fish};
    
    %%
    betas = Betas{i_lr,i_fish};
    b1 = betas(:,1);
    b2 = betas(:,2);
%     b3 = betas(:,3);
    %%
    % set up plot dimensions
    caseflag = 2;
    switch caseflag
        case 1
            X = b1;
            Y = b2;
            Xname = 'motor only (b1)';
            Yname = 'SMT (b2)';
        case 2
            X = b1;%b3;%b1;
            Y = betas(:,5);
%             Y = sqrt(b2.^2+b3.^2);%b2;
            Xname = 'motor only (b1)';
            Yname = 'periodic';
    end
    
    numcell = length(b1);
    
    A = X;
    topN = length(cIX_int);%length(M_cIX{2});%%round(0.01*numcell); % top 5% cutoff
    [~,IX] = sort(A,'descend');
    thresA = A(IX(topN));
    
    B = Y;
    topN = length(cIX_int);%length(M_cIX{2});%round(0.01*numcell); % top 5% cutoff
    [~,IX] = sort(B,'descend');
    thresB = B(IX(topN));%min(Y(cIX_int));%
    
    IX_passX = setdiff(find(A>=thresA),find(B>=thresB));
    IX_passY = setdiff(find(B>=thresB),find(A>=thresA));%find(B>=thresB);%
    IX_pass = union(find(A>=thresA),find(B>=thresB));
    IX_fail = intersect(find(A<thresA),find(B<thresB));%find(A<thresA);
    
    % get min/max
    x0 = min(X(IX_pass));
    x1 = max(X(IX_pass));
    y0 = min(Y(IX_pass));
    y1 = max(Y(IX_pass));
    
    gIX_in = (1:length(X))';
    
    PassX_2{i_lr} = IX_passX;
    PassY_2{i_lr} = IX_passY;
    %% make colormap
    
        clr1 = [0.3,0.7,0.2];
        clr2 = [0.1,0.3,0.9];%[0.3,0.8,1];
        clr_fail = [0.5,0.5,0.5];
%         clr_int_raw = [1,0.2,0.2];%[1,0.1,0];
        clr_int = [0.7,0.2,0.2];
        clrmap = ones(numcell,3);%MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
        clrmap(IX_passX,:) = clr1.*ones(length(IX_passX),3);
        clrmap(IX_passY,:) = clr2.*ones(length(IX_passY),3);
        clrmap(IX_fail,:) = clr_fail.*ones(length(IX_fail),3);
%         clrmap(cIX_int_raw,:) = clr_int_raw.*ones(length(cIX_int_raw),3);
        clrmap(cIX_int,:) = clr_int.*ones(length(cIX_int),3);
             
%     clr1 = [0.3,0.8,0.2];
%     clr2 = [0.1,0.3,1];%[0.3,0.8,1];
%     clr_fail = [0.5,0.5,0.5];
%     clr_int = [1,0.1,0];
%     clrmap = ones(numcell,3);%MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
%     clrmap(IX_passX,:) = clr1.*ones(length(IX_passX),3);
%     clrmap(IX_passY,:) = clr2.*ones(length(IX_passY),3);
%     clrmap(IX_fail,:) = clr_fail.*ones(length(IX_fail),3);
    
    %% bubble plot
    h = figure('Position',[500,100,300,250]); hold on
    U_size = ones(size(X));
    scatter(X,Y,U_size,clrmap,'filled')
    plot([x0,x1],[thresB,thresB],'k--');
    plot([thresA,thresA],[y0,y1],'k--');
    %         plot([x0,x1],[y0,y0],'k--');
    
    scatter(X(IX_passX),Y(IX_passX),1,clr1);%,'filled');
    scatter(X(IX_passY),Y(IX_passY),1,clr2);%,'filled');
    
    scatter(X(cIX_int),Y(cIX_int),1,clr_int_raw);%[1,0.5,0.5]);%,'filled');    
    
%     xlabel(Xname);ylabel(Yname);
    axis equal
    
    xlim([-0.4,0.8]);
    ylim([0,1]);
    set(gca,'XTick',[-0.4,0,0.8],'YTick',[0,1]);
    ylabel('sensory (periodic)');xlabel('motor (m.res.)');
end

%% Plot histogram: all fish pooled, not used in main figure...
isPlotGray = 0;

clr_XY_SEM = [0.7,0.8,0.7; 0.7,0.7,0.8];
% clr_XY = [0.3,0.8,0.2; 0.1,0.3,1];

for i_ax = 1:2
    N_pass = squeeze(N_XY_pass(:,:,i_ax));
    N_int = squeeze(N_XY_int(:,:,i_ax));
    N_full = squeeze(N_XY_full(:,:,i_ax));
    
    %% plot overlapping histogram
    if ~isPlotGray
        figure('Position',[500,500,300,150]);hold on;
    else
        figure('Position',[500,100,100,250]);hold on;
    end
    % plot population average
    meanN = mean(N_pass,1);
    h1 = histogram('BinEdges',edges,'BinCounts',mean(N_pass,1));
    h2 = histogram('BinEdges',edges,'BinCounts',mean(N_int,1));
    
    h1.FaceColor = clr_XY(i_ax,:);
    h1.EdgeColor = 'w';
    h1.FaceAlpha = 0.6;
    h2.FaceColor = [1,0.5,0.5];
    h2.EdgeColor = 'w';
    h2.FaceAlpha = 0.6;
    
    xlim([0,1])
    xlabel('% beta')
    ylabel('a.u.')
    set(gca,'YTick',[]);
    
    % plot SEM in lighter shades
    semN = std(N_pass,0,1)/sqrt(length(range_fish));
    h1_ = histogram('BinEdges',edges,'BinCounts',semN+mean(N_pass,1));
    semN = std(N_int,0,1)/sqrt(length(range_fish));
    h2_ = histogram('BinEdges',edges,'BinCounts',semN+mean(N_int,1));
    
    h1_.FaceColor = clr_XY_SEM(i_ax,:);
    h1_.LineStyle = 'none';
    h1_.FaceAlpha = 0.3;
    h2_.FaceColor = [0.8,0.7,0.7];
    h2_.LineStyle = 'none';
    h2_.FaceAlpha = 0.3;
    
    if isPlotGray
        % plot all cells in gray, safely ignore SEM
        h3 = histogram('BinEdges',edges,'BinCounts',mean(N_full,1));
        
        h3.FaceColor = [0.4,0.4,0.4];
        h3.LineStyle = 'none';
        h3.FaceAlpha = 0.3;
    end
    
end

%% example functional traces
i_fish = i_fish_example;
stimrange = [1,2];
ClusterIDs = [2,1];
cIX_all = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange,1);
%%
% [~,IX]=sort(Y,'descend'); % sort(X)
% i_ex = setdiff(IX(10:20),cIX_int); % intersect

i_stim1_notint = [38927]; % [14940;20564;23555;55905;17715]
i_stim2_int = [20259];%38927; % [13561;16101;16160;18342;18397]
i_int = [12284,15730];%7208; % int_sens [5924;7208;7231;8017;8557;8692]
i_motor = [31789];%7668;% int_motor [6165;6174;6185;6215;7668]
i_motor_notint = [21845];

i_ex = [i_stim2_int,i_int,i_motor];%[i_stim1,i_stim2,i_int,i_motor];
% i_ex = cIX_int(101:100:1000);

h = figure('Position',[500,100,300,250]); hold on
    U_size = ones(size(X));
    scatter(X,Y,U_size,clrmap,'filled')
    plot([x0,x1],[thresB,thresB],'k--');
    plot([thresA,thresA],[y0,y1],'k--');
    %         plot([x0,x1],[y0,y0],'k--');
    
    scatter(X(IX_passX),Y(IX_passX),1,clr1);%,'filled');
    scatter(X(IX_passY),Y(IX_passY),1,clr2);%,'filled');
    
    scatter(X(cIX_int),Y(cIX_int),1,clr_int_raw);%[1,0.5,0.5]);%,'filled');    
    
    clrs = hsv(length(i_ex));
    scatter(X(i_ex),Y(i_ex),50,clrs,'filled');
%     xlabel(Xname);ylabel(Yname);
    axis equal
    
    xlim([-0.4,0.8]);
    ylim([0,1]);
    set(gca,'XTick',[-0.4,0,0.8],'YTick',[0,1]);
    ylabel('sensory (periodic)');xlabel('motor (m.res.)');


    
    % left plot
    cIX = i_ex;
    gIX = (1:length(cIX))';
figure('Position',[50,100,400,500]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',0);
UpdateIndices_Manual(hfig,cIX,gIX);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);


%% multifish bar plot: dividing by anat (midbrain, aHb and pHb)
M_anat_count = zeros(length(range_fish),1);
for i_fishcount = 1:length(range_fish)
    i_fish = range_fish(i_fishcount);
    LoadSingleFishDefault(i_fish,hfig,[1,1],1,0);
    %% divide by anat
    cIX_int = Intersect_fixedprct{i_fish};
    
    cIX = cIX_int;
    gIX = ones(size(cIX));
    
    % count cells: hindbrain vs not
    MASKs = getappdata(hfig,'MASKs');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    
    %     Msk_IDs = [94,219,220];% midbrain 94; Rh1 219; Rh2 220; hindbrain 114;
    cIX_mb = ScreenCellsWithMasks(94,cIX,gIX,MASKs,CellXYZ_norm,absIX);
    cIX_ahb = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
    cIX_phb = ScreenCellsWithMasks([221:225],cIX,gIX,MASKs,CellXYZ_norm,absIX);
    
    M_anat_count(i_fishcount,1) = length(cIX_mb)/length(absIX);
    M_anat_count(i_fishcount,2) = length(cIX_ahb)/length(absIX);
    M_anat_count(i_fishcount,3) = length(cIX_phb)/length(absIX);
    
end
%% multi-fish bar plot of the above
figure('Position',[100,400,150,160]);hold on;

inc1 = 0.2;
inc = 0.15;

% connect the x columns with grey lines for each fish
for i_fishcount = 1:length(range_fish)
    Y = M_anat_count(i_fishcount,:)*100;
    plot([1,2,3],Y,'color',[0.7,0.7,0.7],'linewidth',0.5)
end

% plot the data points as dots, with avr and SEM
for x = 1:3
    Y = M_anat_count(:,x)*100;
    
    
    %     scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
    scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
        0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
    plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
    err = std(Y)/sqrt(length(Y));
    plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
    plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
    plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);
end

xlim([0.5,3.5])
set(gca,'XTick',[1,2,3],'XTickLabels',{'midbrain','hindbrain Rh1,2','hindbrain Rh3+'},'XTickLabelRotation',45);
ylabel('% of all cells')

%% bar plots for multifish data (obsolete)
% to show that difference between int & motor is significant
% but actually difference between int & sensory is significant too just not
% as big
for i_ax = 1:2
    
    figure('Position',[100,400,150,160]);hold on;
    
    inc1 = 0.2;
    inc = 0.15;
    
    for x = 1:3
        switch x
            case 1
                Y = squeeze(s_XY_pass(:,1,i_ax));
            case 2
                Y = squeeze(s_XY_int(:,1,i_ax));
            case 3
                Y = squeeze(s_XY_full(:,1,i_ax));
        end
        scatter(x*ones(length(Y),1),Y,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',...
            0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
        err = std(Y)/sqrt(length(Y));
        plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
        plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);
        plot([x,x],[mean(Y)-err,mean(Y)+err],'color',[1,0.5,0.5]);
    end
    
    xlim([0.5,3.5])
    if i_ax==1
        set(gca,'XTickLabels',{'motor','convergent','all'},'XTickLabelRotation',45);
    else
        set(gca,'XTickLabels',{'sensory','convergent','all'},'XTickLabelRotation',45);
    end
end
% ylabel('')

%% t-test
i_ax = 1; % motor
y1 = squeeze(s_XY_pass(:,1,i_ax));
y2 = squeeze(s_XY_int(:,1,i_ax));
[h,p] = ttest2(y1,y2)

i_ax = 2; % sensory
y1 = squeeze(s_XY_pass(:,1,i_ax));
y2 = squeeze(s_XY_int(:,1,i_ax));
[h,p] = ttest2(y1,y2)
