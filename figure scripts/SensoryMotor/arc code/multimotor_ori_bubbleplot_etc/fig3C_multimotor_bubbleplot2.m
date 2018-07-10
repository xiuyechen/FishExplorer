%% Autoclus based analysis
% Load data in GUI, rank Autoclusters by multimotor, and export to workspace




%% 
% %% threshold the bubble plot with chosen radius
% thres_rad = 0.3;
% 
% radius = sqrt(stimcorr.^2 + motorcorr.^2);
% [A,U_sorted] = sort(radius,'descend');
% i_end = find(A>thres_rad,1,'last');
% disp(i_end)
% IX_passrad = U_sorted(1:i_end);
% 
% cIX_radthres = []; gIX_radthres = [];
% for i = 1:length(IX_passrad),
%     IX = find(gIX_in==IX_passrad(i));
%     cIX_radthres = [cIX_radthres; cIX_in(IX)];
%     gIX_radthres = [gIX_radthres; gIX_in(IX)];
% end
% %% Anat plot of only clusters > radius, same colormap
% isRefAnat = 1;
% isPopout = 1;
% figure
% DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_radthres,gIX_radthres,clrmap);
% 
% %% bubble plot in 2-D color, with radius drawn
% figure('Position',[500,500,300,250]);hold on;
% scatter(stimcorr,motorcorr,U_size,clrmap)
% xlabel('stimulus corr.');ylabel('motor corr.');
% theta = -1:0.01:pi/2;
% X = cos(theta)*thres_rad;
% Y = sin(theta)*thres_rad;
% plot(X,Y,'k--','Linewidth',1.5)
% axis equal
% % ylim([-0.15,0.4])
% % xlim([0,0.6])
% xlim([0,0.7]);
% ylim([-0.22,0.6]);
% set(gca,'YTick',-0.2:0.2:0.6);
% 
% %%
% thres_x = 0.3;% 0.1; for fish8 figure
% [A,U_sorted] = sort(stimcorr,'ascend');
% i_end = find(A<thres_x,1,'last');
% IX_smallstim = U_sorted(1:i_end);
% IX_plot = intersect(IX_smallstim,IX_passrad,'stable');
% disp(length(IX_plot));
% 
% cIX_out = []; gIX_out = [];
% for i = 1:length(IX_plot),
%     IX = find(gIX_in==IX_plot(i));
%     cIX_out = [cIX_out; cIX_in(IX)];
%     gIX_out = [gIX_out; i*ones(size(IX))];
% end
% %%
% cIX = cIX_out;
% gIX = gIX_out;

%% multi-motor on individual cells
i_fish = 8;
isFullData = true;
LoadFullFish(hfig,i_fish,isFullData);

setappdata(hfig,'isZscore',1);
    
i_ClusGroup = 6;
i_Cluster = 1;

[cIX_load,gIX_load] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
% cIX = cIX_load(1:10:end);
% gIX = (1:length(cIX))';%gIX_load(1:10:end);

M_stimrange = GetStimRange(); % get the default stim-range (no input param)   
stimrange = M_stimrange{i_fish};
timelists = getappdata(hfig,'timelists');
fishset = getappdata(hfig,'fishset');

tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
[M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);

cIX = cIX_load;
gIX = gIX_load;

%% 
% first, export from workspace
gIXCopy = gIX;

%% cell based
gIX = (1:length(cIX))';
tic
[betas,stimcorr,motorcorr] = MultiMotorRegression(i_fish,M,stim,behavior);
toc
MultiMotorVisuals(hfig,betas,cIX,gIX);

%% cluster based
gIX = gIXCopy;
C = FindClustermeans(gIX,M);
tic
[betas_C,stimcorr,motorcorr] = MultiMotorRegression(i_fish,C,stim,behavior);
toc
MultiMotorVisuals(hfig,betas_C,cIX,gIX);

%%
% figure;subplot(121);hist(betas_C(:,end),-20:0.5:40);subplot(122);hist(betas(:,end),-20:0.5:40);

%% threshold and visualize
% thres_stim = prctile(stimcorr,90);
% thres_motor = prctile(motorcorr,90);
% IX_stim = find(stimcorr>thres_stim);
% IX_motor = find(motorcorr>thres_motor);
% IX = union(IX_stim,IX_motor);
% MultiMotorVisuals(hfig,betas(IX,:),cIX(IX),gIX(IX));



