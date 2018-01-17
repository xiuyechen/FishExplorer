%%
clear all; close all; clc

%% load
outputDir = GetOutputDataDir;
load(fullfile(outputDir,'PTintOMR_3%each.mat'),'PTintOMR');
load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));

range_fish = 8:18;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
XY_LR = cell(2,2,18);
XY_LR_IXpass = cell(2,2,18);
for i_fish = range_fish
    cIX_int = PTintOMR{i_fish};
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
                A = sqrt(b2.^2+b3.^2);%b2;;
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

for i_fishcount = 1:length(range_fish)
    i_fish = range_fish(i_fishcount);
    cIX_int = PTintOMR{i_fish};
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
            if true%false
                % i_fish=8, i_lr=2;
                figure('Position',[500,500,300,150]);hold on;
                
                h1 = histogram(y_pass_n);
                h2 = histogram(y_int_n);
                
                h1.Normalization = 'probability';
                h1.BinWidth = 0.025;
                h2.Normalization = 'probability';
                h2.BinWidth = 0.025;
                
                h1.FaceColor = clr_XY(i_ax,:);
                h1.EdgeColor = clr_XY(i_ax,:);
                h2.FaceColor = [1,0.5,0.5];
                h2.EdgeColor = [1,0.5,0.5];
                
                xlim([0,1])
                xlabel('% beta')
                ylabel('a.u.')
                %     ylim([0,inf])
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

%% bar plots for multifish data (obsolete)
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
