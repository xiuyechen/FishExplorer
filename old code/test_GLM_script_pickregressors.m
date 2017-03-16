% temp script for GLM, fig5

isSaveFig = 1;
isPlotFig = 1;

outputDir = GetOutputDataDir;
saveDir = [];
saveDir{1} = fullfile(outputDir,'fig5_GLM_ranked_Left_030917');
saveDir{2} = fullfile(outputDir,'fig5_GLM_ranked_Right_030917');
if ~exist(saveDir{1}, 'dir'), mkdir(saveDir{1}), end;
if ~exist(saveDir{2}, 'dir'), mkdir(saveDir{2}), end;

%%
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%

thres_excl = 0.6;

range_fish = GetFishRange;
for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    
%     setappdata(hfig,'isTrialRes',1);
    
    
    %% load all Autoclus centroids as full regressor set
    ClusterIDs = [6,1];
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    C = FindClustermeans(gIX_load,M);
    Cpos = ShiftBaselinePosSignal(C',1);
    
    nClus = size(Cpos,2);
    Data = behavior;
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
    for i_lr = 1%:2
        y = Data(i_lr,:)';
        
        Rcorr_motor = corr(Cpos,y);      
        
        %%
        rIX_keep = [];
        for i_itr = 1:nClus % dummy index
            disp(i_itr);
            rIX_test0 = setdiff(1:nClus,rIX_keep);
            
            M_r2 = nan(nClus,1);
            
            rIX_tested = [];
            for i_clus = rIX_test0
                % exclude 0.7 correlated clusters
                X_old = Cpos(:,rIX_keep);
                X_avoid = [X_old,y];
                X_test = Cpos(:,i_clus);
                R = corr(X_avoid,X_test);
                
                if max(R)<thres_excl
                    rIX_tested = [rIX_tested,i_clus]; %#ok<AGROW>
                    X_wConst = [ones(size(Cpos,1),1),X_old,X_test];
                    [~,~,~,~,stat] = regress(y,X_wConst);
                    M_r2(i_clus,1) = stat(1);
                end
            end
            
            if ~isempty(rIX_tested)
                % end of round, pick best regressor to add to collection
                H = abs(Rcorr_motor(rIX_tested));
                [~,ix] = min(H);
                rIX_new = rIX_tested(ix);

%                 % calculated based on 2 factors: var explained, and (low) motor score
%                 x_VE = M_R2(rIX_tested);
%                 y_H = 1-abs(Rcorr_motor(rIX_tested));
%                 Rad = sqrt(x_VE.^2+y_H.^2);
%                 [~,ix] = max(Rad);
%                 rIX_new = rIX_tested(ix);
                
%                 % calculated based on 2 factors: var explained, and stim-lock score
%                 x_VE = M_R2(rIX_tested);
%                 y_H = score_StimLock(rIX_tested);
%                 Rad = sqrt(x_VE.^2+y_H.^2);
%                 [~,ix] = max(Rad);
%                 rIX_new = rIX_tested(ix);
                
                rIX_keep = [rIX_keep,rIX_new];
                M_r2_keep(i_itr) = M_r2(rIX_new);
            end
        end
    
    
    %% scatter-plot: var-expl as function of motor/stimlock
    IX = rIX_keep;
    cmap = flipud(jet(length(IX)));
    h1 = figure('Position',[200,200,600,300]);
    subplot(1,2,1); hold on
    scatter(M_r2_keep(1:length(IX))',abs(Rcorr_motor(IX)),C_size_score(IX),cmap)
    plot([0,1],[thres_excl,thres_excl],'r:');
    maxR2 = M_r2_keep(length(IX));
    plot([maxR2,maxR2],[0,1],'r--');
    text(maxR2-0.4,0.1,['var.expl=',num2str(maxR2,2)]);
    xlabel('var explained');
    ylabel('abs(motor.corr)');
    xlim([0,1]);
    ylim([0,1]);
%     axis equal
    
    subplot(1,2,2); hold on
    [~,~,H] = GetTrialAvrLongTrace(hfig,Cpos');
    VE = M_r2_keep(1:length(IX));
    scatter(VE,score_StimLock(IX),C_size_score(IX),cmap)
    plot([0,1],[1,1],'r:');
    xlabel('var explained');
    ylabel('stimlock score');
    xlim([0,1]);
    ylim([0,1]);
%     axis equal
    
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
    gIX_rank = rIX_keep;
    %%
    % set rank
    [gIX,cIX] = SetGroupIXorder(gIX_load,gIX_rank,cIX_load);
    UpdateIndices_Manual(hfig,cIX,gIX);
    
    h2 = figure('Position',[50,100,1400,800]);
    % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
    subplot(121)
    %     setappdata(hfig,'isPlotBehavior',1);
    %     setappdata(hfig,'isStimAvr',0);
    setappdata(hfig,'isPlotLines',0);
    setappdata(hfig,'clrmap_name','jet');
    %     UpdateTimeIndex(hfig);
    
    DrawTimeSeries(hfig,cIX,gIX);
    
    % right plot
    ax = subplot(122)
    I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
    DrawCellsOnAnat(I,ax);
    
    %% Save plot
    if isSaveFig
        f = combineFiguresLR(h1, h2);
        filename = fullfile(saveDir{i_lr}, ['Fish',num2str(i_fish)]);
        saveas(gcf, filename, 'png');
        close(gcf)
    end
    end
end

%% view stim/motor in 2D scatter plot
figure;
cmap = hsv(nClus);
scatter(score_StimLock,Rcorr_motor,C_size_score,cmap);

%% motor: tAvr/total variance
[bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior);
vAvr = var(bh_tAvr(1,:));
vTot = var(behavior(1,:));

vRes = var(bh_tRes(1,:));

%% Leave-1-out of current model: not very effective
X_all = Cpos(:,rIX_keep);
nReg = length(rIX_keep);

M_KO_r2 = zeros(nReg,1);
for i_reg = 1:nReg
    IX = find((1:nReg)~=i_reg);
    X = X_all(:,IX);
    X_wConst = [ones(size(X,1),1),X];
    [~,~,~,~,stat] = regress(y,X_wConst);
    M_KO_r2(i_reg,1) = stat(1);
end

figure;
plot(sort(M_KO_r2));

%% Regularization
tic
X = Cpos(:,rIX_keep);
X_wConst = [X,ones(size(X,1),1)];
[B0,FitInfo] = lasso(X_wConst,y,'CV',2);
toc

%%

lassoPlot(B0,FitInfo,'PlotType','CV');
B = B0(2:end,51);
%%
[B_sorted,IX] = sort(B,'descend');
rIX_ordered = rIX_keep(IX);
[gIX,cIX] = SetGroupIXorder(gIX_load,rIX_ordered,cIX_load);

%%
numC = 64;
cmap = Make1DColormap([0,0,1;1,1,1;1,0,0],numC);
Xrange = [-0.1,0.1];
clrIX = MapXto1Dcolormap(B_sorted,Xrange,numC);
clrmap = cmap(clrIX,:);

UpdateIndices_Manual(hfig,cIX,gIX);
    
    h2 = figure('Position',[50,100,1400,800]);
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
% B = lasso(X,y,'Lambda',vLambda);% 'CV',2);
%%
% x_VE = M_R2;
%
% [~,~,H] = GetTrialAvrLongTrace(hfig,C);
% y_H = 1-H;
% Rad = sqrt(x_VE.^2+y_H.^2);
%
% %%
% IX = find(M_R2>0.2);
% figure;
% C_plot = zscore(C(IX,:),0,2);
% bh = zscore(behavior,0,2);
% imagesc(vertcat(C_plot,bh))
%
% corr(y,C(IX,:)')