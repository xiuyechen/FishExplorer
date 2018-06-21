% temp script for GLM, fig5
clear all;close all;clc
%%
isSaveFig = 1;

outputDir = GetOutputDataDir;
saveDir = [];
saveDir{1} = fullfile(outputDir,'fig5_GLM_L_031417');
saveDir{2} = fullfile(outputDir,'fig5_GLM_R_031417');
if ~exist(saveDir{1}, 'dir'), mkdir(saveDir{1}), end;
if ~exist(saveDir{2}, 'dir'), mkdir(saveDir{2}), end;

%%
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
thres_excl = 0.5;

range_fish = GetFishRange;
for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);

    %% load all Autoclus centroids as full regressor set
    ClusterIDs = [6,1];
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    C = FindClustermeans(gIX_load,M);
    Cpos = ShiftBaselinePosSignal(C',1);
    
    nClus = size(Cpos,2);
    Data = behavior';
    [~,~,H] = GetTrialAvrLongTrace(hfig,Cpos');
    score_StimLock = 1-H;
    
    U = unique(gIX_load);
    C_size = zeros(nClus,1);
    for i=1:length(U)
        C_size(i) = length(find(gIX_load==U(i)));
    end
    C_size_score = sqrt(C_size);
    
    %%
    
    M_r2_keep = nan(nClus,1);
    for i_lr = 1:2
        %% get behavior trace to model
        y = Data(:,i_lr);
        
        %% exclude clusters that are too similar to behavior
        Rcorr_motor = max(corr(Cpos,Data),[],2);
        [R_sort,IX] = sort(abs(Rcorr_motor));
        ix_end = find(R_sort<thres_excl,1,'last');
        rIX_cand = IX(1:ix_end);
        %% step-wise include more clusters into model, get r2(var.expl)
        rIX_select = [];
        for i_itr = 1:ix_end % dummy index
            rIX_this = rIX_cand(i_itr);
            X_add = Cpos(:,rIX_this);
            X_old = Cpos(:,rIX_select);
            X_wConst = [ones(size(Cpos,1),1),X_old,X_add];
            
            [~,~,~,~,stat] = regress(y,X_wConst);
            
            rIX_select = [rIX_select,rIX_this];
            M_r2_keep(i_itr) = stat(1);
        end

        %% scatter-plot: var.expl as function of motor (/stimlock see end)
        IX = rIX_select;
        cmap = flipud(jet(length(IX)));
        f = [];
        f{1} = figure('Position',[200,200,300,300]);
        hold on
        scatter(abs(Rcorr_motor(IX)),M_r2_keep(1:length(IX))',C_size_score(IX),cmap)
        plot([thres_excl,thres_excl],[0,1],'r:');
        maxR2 = M_r2_keep(length(IX));
        plot([0,1],[maxR2,maxR2],'r--');
        text(0.55,maxR2+0.05,['var.expl=',num2str(maxR2,2)]);
        ylabel('var explained');
        xlabel('abs(motor.corr)');
        ylim([0,1]);
        xlim([0,1]);

        %% Lasso regularization to eliminate mostly redundant clusters
        tic
        X = Cpos(:,rIX_select);
        X_wConst = [ones(size(X,1),1),X];
        [B0,FitInfo] = lasso(X,y);%,'CV',10);
        toc

        %% get coeff of determination 'r2'
        r2s = 1 - FitInfo.MSE/var(y);
        
        %% select model with a particular regularization strength/lambda 
        % current criterion: only leave ~15 clusters
        % by default, the lasso code tests 100 lambda values
        DF = FitInfo.DF;
        ix = find(DF<=15,1,'first');
        % ix = FitInfo.Index1SE; % e.g. for default lasso with CV; conservative 

        % (code below: for this example lambda only)
            
        %% show trends with strength of regularization
        f{2} = figure('Position',[300,300,300,200]);
        hold on;
        yyaxis left
        plot(r2s)
        ylabel('r2')
        plot([ix,ix],[0,1],'r--')
        text(ix+1,r2s(ix)+0.05,num2str(r2s(ix),2))
        
        yyaxis right
        plot(FitInfo.DF)
        ylabel('#of dims')
        % ylim([-150 150])
        grid on
        xlabel('strength of regularization')

        %% compare modeled behavior with original 
        f{3} = figure('Position',[200,200,700,200]);
        hold on;
        plot(y+8,'r')
        plot(X*B0(:,ix)+FitInfo.Intercept(ix));
        text(-300,8,'orig.')
        text(-300,0,'model')
        axis off
  
        %%
        B = B0(1:end,ix);
        [B_sorted,IX] = sort(B,'descend');
        rIX_sorted = rIX_select(IX);
        
        thres_plot = prctile(abs(B_sorted),80);
        ix1 = find(B_sorted>thres_plot);
        ix2 = find(B_sorted<-thres_plot);
        IX = [ix1;ix2];
        rIX_plot = rIX_sorted(IX);
        B_plot = B_sorted(IX);
        [gIX,cIX] = SetGroupIXorder(gIX_load,rIX_plot,cIX_load);
        
        %%
        numC = 64;
        cmap = Make1DColormap([0,0,1;1,1,1;1,0,0],numC);
        Xrange = [-0.2,0.2];
        clrIX = MapXto1Dcolormap(B_plot,Xrange,numC);
        clrmap = cmap(clrIX,:);
        
        UpdateIndices_Manual(hfig,cIX,gIX);
        
        f{4} = figure('Position',[50,100,1400,800]);
        % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
        subplot(121);
        %     setappdata(hfig,'isPlotBehavior',1);
        %     setappdata(hfig,'isStimAvr',0);
        setappdata(hfig,'isPlotLines',0);
        %     setappdata(hfig,'clrmap_name','jet');
        %     UpdateTimeIndex(hfig);
        opts = []; opts.clrmap = clrmap;
        DrawTimeSeries(hfig,cIX,gIX,opts);
        
        % right plot
        ax = subplot(122);
        I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
        DrawCellsOnAnat(I,ax);

        %% get top PC's from this set of (clusmean) regressors
        % do pca
        [coeff,score,latent,tsquared,explained,mu] = pca(X'); % ~52 sec on linux
        
        % plot top PC's, weighted by var.explained
        f{5} = figure('Position',[300,300,300,200]);
        n_top = 30;                
        im = coeff;
        im = im.*repmat(explained',size(coeff,1),1);
        imagesc(im(:,1:n_top))
        xlabel('PC''s, weighted')
       
        topPCs = coeff(:,1:n_top);
        
        % fit model (for motor) with increasing number of PC's, and calculate r2 stat
        nPCs = size(coeff,2);
        r2s_numPC = zeros(nPCs,1);
        for i = 1:nPCs
            topPCs = coeff(:,1:i);
            X_pc = [ones(size(topPCs,1),1),topPCs];
            [b,~,~,~,stat] = regress(y,X_pc);
            r2s_numPC(i) = stat(1);
        end
        
        % plot r2 as function of number of PC's included in model
        f{6} = figure('Position',[300,300,300,200]);
        plot(r2s_numPC)
        grid on
        ylabel('var explained (r2)')
        xlabel('# of PC''s')
        
        %% Save plot
        h1 = combineFiguresLR(f{1},f{2});
        h1 = combineFiguresTB(h1,f{3});
        h2 = combineFiguresLR(f{5},f{6});        
        h = combineFiguresTB(h1,h2);
        
        SaveFigureHelper(isSaveFig, saveDir{i_lr}, ['Fish',num2str(i_fish)]);
        SaveFigureHelper(isSaveFig, saveDir{i_lr}, ['Fish',num2str(i_fish),'_anat'],f{4});

    end
end


        %% 'CV'
        
        % lassoPlot(B0,FitInfo,'PlotType','CV');
        % B = B0(2:end,ix);
        %
        % FitInfo.MSE(FitInfo.Index1SE)
          
        %% r2
       %         r2s = zeros(100,1);
%         for i = 1:100;
%             y_hat = X*B0(:,i)+FitInfo.Intercept(i);
%             SSE = sum((y-y_hat).^2);
%             SST = sum((y-mean(y)).^2);
%             r2s(i) = 1-SSE/SST;
%         end
%         FitInfo.r2 = r2s;
               
  %%
        %     dVE = diff([0;VE]);
        %     figure;hold on
        %     scatter(dVE,score_StimLock(IX),C_size_score(IX),cmap)
        % %     plot([0,1],[1,1],'r:');
        %     xlabel('var explained');
        %     ylabel('stimlock score');
        % %     xlim([0,1]);
        % %     ylim([0,1]);
        % %     axis equal
        % [dVE_sort,IX_sort] = sort(dVE,'descend');
        % ix_end = 20;%find(dVE_sort>0.02,1,'last');
        %     gIX_rank = IX_sort(1:ix_end);
        
        %%
        %         gIX_rank = rIX_select;
        %         %%
        %         % set rank
        %         [gIX,cIX] = SetGroupIXorder(gIX_load,gIX_rank,cIX_load);
        %         UpdateIndices_Manual(hfig,cIX,gIX);
        %
        %         h2 = figure('Position',[50,100,1400,800]);
        %         % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
        %         subplot(121)
        %         %     setappdata(hfig,'isPlotBehavior',1);
        %         %     setappdata(hfig,'isStimAvr',0);
        %         setappdata(hfig,'isPlotLines',0);
        %         setappdata(hfig,'clrmap_name','jet');
        %         %     UpdateTimeIndex(hfig);
        %
        %         DrawTimeSeries(hfig,cIX,gIX);
        %
        %         % right plot
        %         ax = subplot(122);
        %         I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
        %         DrawCellsOnAnat(I,ax);
        %
        
        
        
        %% view stim/motor in 2D scatter plot
        % figure;
        % cmap = hsv(nClus);
        % scatter(score_StimLock,Rcorr_motor,C_size_score,cmap);
        
        
        %% Leave-1-out of current model: not very effective
        % X_all = Cpos(:,rIX_select);
        % nReg = length(rIX_select);
        %
        % M_KO_r2 = zeros(nReg,1);
        % for i_reg = 1:nReg
        %     IX = find((1:nReg)~=i_reg);
        %     X = X_all(:,IX);
        %     X_wConst = [ones(size(X,1),1),X];
        %     [~,~,~,~,stat] = regress(y,X_wConst);
        %     M_KO_r2(i_reg,1) = stat(1);
        % end
        %
        % figure;
        % plot(sort(M_KO_r2));
        
        %% plot stim-lock score against var-explained, for reg selection
%          subplot(1,2,2); hold on
%         [~,~,H] = GetTrialAvrLongTrace(hfig,Cpos');
%         VE = M_r2_keep(1:length(IX));
%         scatter(VE,score_StimLock(IX),C_size_score(IX),cmap)
%         plot([0,1],[1,1],'r:');
%         xlabel('var explained');
%         ylabel('stimlock score');
%         xlim([0,1]);
%         ylim([0,1]);
        