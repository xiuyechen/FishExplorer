%% motor: tAvr/total variance
behavior_z = zscore(behavior,0,2);
[bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior_z);
vAvr_bh = var(bh_tAvr(1,:));
vTot_bh = var(behavior_z(1,:));

vRes_bh = var(bh_tRes(1,:));

%% make histogram for all cells
M_ = M_0(1:100:end,:);
[M_tAvr,M_tRes] = GetTrialAvrLongTrace(hfig,M_);
vAvr = var(M_tAvr,0,2);
vRes = var(M_tRes,0,2);
vTot = var(M_,0,2);
figure;
histogram(vAvr)
% figure;
% histogram(vRes)

% Rcorr = max(corr(M_',behavior'),[],2);
%%
% behavior = 
Rcorr = corr(M_',bh_tRes(1,:)');

%%
figure('Position',[300,300,200,200]); hold on
scatter(vRes,Rcorr,ones(nCells,1)*0.1);
axis equal
xlim([0,1]);
ylim([0,1]);
xlabel('var.tRes');
ylabel('corr to motorseed')
title('left motorseed')
%%
nCells = size(M_0,1);
clr = repmat([0.7,0.7,0.7],nCells,1);
absIX = getappdata(hfig,'absIX');
nClus = length(unique(gIX_load));
clrmap = hsv(nClus+1);
clr(cIX_load,:) = clrmap(gIX_load,:);
%
figure; hold on
scatter(vRes,Rcorr,ones(nCells,1),clr);
axis equal
xlim([0,1]);
ylim([0,1]);
xlabel('var.tRes');
ylabel('corr to motorseed')
title('left motorseed')
%%
nClus = length(unique(gIX_load));
clrmap = hsv(nClus+1);
clr = clrmap(gIX_load,:);
%
figure; hold on
scatter(vAvr(cIX_load),Rcorr(cIX_load),ones(length(cIX_load),1),clr);
axis equal
xlim([0,1]);
ylim([0,1]);
xlabel('var.tAvr');
ylabel('corr to R-motorseed')
title('left motorseed')