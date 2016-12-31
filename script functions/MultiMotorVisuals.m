function MultiMotorVisuals(hfig,betas,cIX_in,gIX_in)
%% get coefficients from multiple linear regression
% cIX_in = getappdata(hfig,'cIX');
% gIX_in = getappdata(hfig,'gIX_betas');
% betas = getappdata(hfig,'betas');
stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
% figure;scatter(stimcorr,motorcorr)%scatter(motorcorr,stimcorr)

%% draw custom 2-D colormap illustration (square)
% res = 100;
% grad = linspace(0,1,res);
% rev_grad = linspace(1,0,res);
% 
% grid = ones(res,res,3);
% 
% grid(:,:,3) = repmat(grad,res,1)';
% grid(:,:,1) = 0.5*repmat(rev_grad',1,res)'+0.5*repmat(rev_grad,res,1)';
% grid(:,:,2) = repmat(grad',1,res)';
% 
% clrmap_2D = reshape(grid,res*res,3);
% 
% % figure;imagesc(grid)
% % axis xy
% % axis off
% % axis equal
[grid, res] = MakeDiagonal2Dcolormap;

% [plot scaled 2D-colorbar?]

% figure;imagesc(grid)
% axis xy
% axis equal


%% get new gIX with matching custom colormap 'cmap_U'
if false, % temp, trying this for cell-based betas
    thres_stim = prctile(stimcorr,90);
    thres_motor = prctile(motorcorr,90);
    IX_stim = find(stimcorr>thres_stim);
    IX_motor = find(motorcorr>thres_motor);
    stimcorr(IX_stim) = thres_stim;
    motorcorr(IX_motor) = thres_motor;
end

gIX_x = round((stimcorr-min(stimcorr))/(max(stimcorr)-min(stimcorr))*(res-1))+1;
gIX_y = round((motorcorr-min(motorcorr))/(max(motorcorr)-min(motorcorr))*(res-1))+1;

gIX2 = SqueezeGroupIX(gIX_in);
U = unique(gIX2);
U_size = zeros(size(gIX_x));
clrmap = zeros(length(gIX_x),3);
clrmap_2D = reshape(grid,res*res,3); % for efficient indexing

for i = 1:length(U);
    ix = find(gIX2 == U(i));
    U_size(i) = length(ix);
    ix = sub2ind([res,res],gIX_y(U(i))',gIX_x(U(i))');
    clrmap(i,:) = clrmap_2D(ix,:);
end

%% bubble plot in 2-D color (plot of all clusters, cluster size indicated by circular marker size) 
figure('Position',[500,500,300,250]);
scatter(stimcorr,motorcorr,U_size,clrmap)
xlabel('stimulus corr.');ylabel('motor corr.');
axis equal
% xlim([0,0.7]);
% ylim([-0.22,0.6]);
% set(gca,'YTick',-0.2:0.2:0.6);

%% Anat plot with custom colormap
% isRefAnat = 1;
% isPopout = 1;
figure
% DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_in,gIX_in,clrmap);
% DrawCellsOnAnatProj_othercolor(hfig,cIX_in,gIX_in,cmap_U,isRefAnat,isPopout);
opts = [];
opts.isShowFishOutline = true;
opts.isPopout = true;
I = LoadCurrentFishForAnatPlot(hfig,cIX_in,gIX_in,clrmap,[],opts);
DrawCellsOnAnat(I);