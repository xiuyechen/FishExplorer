% fig5: sensory-motor map. 
% fig5A: PT-L vs PT-R, in 2-D colors

% using Fish8 as example for now

% export data for stimrange=1 (PT) only

% Data = getappdata(hfig,'M_0');
% %         Data = getappdata(hfig,'M');
% thres_reg = getappdata(hfig,'thres_reg');
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
% behavior = getappdata(hfig,'behavior');
% % cIX_in = getappdata(hfig,'cIX');
% % gIX_in = getappdata(hfig,'gIX');
% % numK = getappdata(hfig,'numK');
i_fish = getappdata(hfig,'i_fish');

%%
% get stim/motor regressors
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
%%
cIX_plot = cIX_L;
gIX_plot = gIX_L;
clrmap = jet(64);

% [cIX_plot,ia,ib] = union(cIX_L,cIX_R);
% corr_plot_L = [Corr_L(IX_L(ia)),Corr_L(IX_R(ib))];
% corr_plot_R = [Corr_R(IX_L(ia)),Corr_R(IX_R(ib))];

%%

clr1 = [1,1,1];
clr2 = [1,0,0];
numC = 64;
clrmap = makeColormap(clr1,clr2,numC);

cmap = clrmap;%(unique(gIX_plot),:);
figure; 
I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,cmap);
DrawCellsOnAnat(I);

%% interpolate 1-D colormap
m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
bottom = [1 0 0];middle = [0 0 0]; top = [0 1 0];
new = [bottom; middle; top];
% x = 1:m;

oldsteps = linspace(0, 1, 3);
newsteps = linspace(0, 1, m);
newmap = zeros(m, 3);

for i=1:3
    % Interpolate over RGB spaces of colormap
    newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1); % confine between range [0,1]
end

figure
set(gcf, 'colormap', newmap), colorbar

%%


% % red-white-blue colormap
% cmap = zeros(64,3);
% cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
% cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
% cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];