% fig5: sensory-motor map. 
% fig5A: PT-L vs PT-R, in 2-D colors

% using Fish8 as example for now

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 8;
ClusterIDs = GetClusterIDs('all');
stimrange = 1;
[cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);

%% regression for phototaxis
% ResetDisplayParams(hfig,i_fish);

% get stim/motor regressors
stim = getappdata(hfig,'stim');
fishset = getappdata(hfig,'fishset');
[~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);

reg_range = [2,3];
reg_thres = 0.3;

% Left
i_reg = 2;
Reg = regressor_s(i_reg,:);
Corr_L = corr(Reg',M');

IX_L = find(Corr_L>reg_thres);
corr_L = Corr_L(IX_L);
cIX_L = cIX(IX_L);
gIX_L = (ceil((corr_L-0.3)*63/(1.0-0.3))+1)';

% Right 
i_reg = 3;
Reg = regressor_s(i_reg,:);
Corr_R = corr(Reg',M');

IX_R = find(Corr_R>reg_thres);
corr_R = Corr_R(IX_R);
cIX_R = cIX(IX_R);
gIX_R = ceil((corr_R-0.3)*63/(1.0-0.3))+1;

%% Plot left side on red-white scale
cIX_plot = cIX_L;
gIX_plot = gIX_L;

[cIX_plot,ia,ib] = union(cIX_L,cIX_R);
corr_plot_L = [Corr_L(IX_L(ia)),Corr_L(IX_R(ib))];
corr_plot_R = [Corr_R(IX_L(ia)),Corr_R(IX_R(ib))];

%%

clr1 = [1,1,1];
clr2 = [1,0,0];
numC = 64;
clrmap = MakeColormap(clr1,clr2,numC);

cmap = clrmap;%(unique(gIX_plot),:);
figure; 
I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,cmap);
DrawCellsOnAnat(I);

%%


% % red-white-blue colormap
% cmap = zeros(64,3);
% cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
% cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
% cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
