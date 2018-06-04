function [cIX,gIX] = BubblePlot(hfig)
isplotfig = 0;

cIX_in = getappdata(hfig,'cIX');
gIX_in = getappdata(hfig,'gIX_betas');

%% get coefficients from multiple linear regression
betas = getappdata(hfig,'betas');
stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
% figure;scatter(stimcorr,motorcorr)%scatter(motorcorr,stimcorr)

%% draw custom 2-D colormap illustration (square)
res = 100;
grad = linspace(0,1,res);
rev_grad = linspace(1,0,res);

grid = ones(res,res,3);

grid(:,:,3) = repmat(grad,res,1)';
grid(:,:,1) = 0.5*repmat(rev_grad',1,res)'+0.5*repmat(rev_grad,res,1)';
grid(:,:,2) = repmat(grad',1,res)';

clrmap_2D = reshape(grid,res*res,3);

% figure;imagesc(grid)
% axis xy
% axis off
% axis equal

%% get new gIX with matching custom colormap 'cmap_U'
gIX_x = round((stimcorr-min(stimcorr))/(max(stimcorr)-min(stimcorr))*(res-1))+1;
gIX_y = round((motorcorr-min(motorcorr))/(max(motorcorr)-min(motorcorr))*(res-1))+1;

gIX_old = gIX_in;
U = unique(gIX_old);
U_size = zeros(size(gIX_x));
clrmap = zeros(length(gIX_x),3);
for i = 1:length(U);
    ix = find(gIX_old == U(i));
    U_size(i) = length(ix);
    ix = sub2ind([res,res],gIX_y(U(i))',gIX_x(U(i))');
    clrmap(i,:) = clrmap_2D(ix,:);
end

%% bubble plot in 2-D color (plot of all clusters, cluster size indicated by circular marker size)
% figure('Position',[500,500,250,200]);scatter(stimcorr,motorcorr,U_size,clrmap)
% xlabel('stimulus corr.');ylabel('motor corr.');

%% Anat plot with custom colormap
if isplotfig,
    isRefAnat = 1;
    isPopout = 1;
    figure
    DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_in,gIX_in,clrmap);
    % DrawCellsOnAnatProj_othercolor(hfig,cIX_in,gIX_in,cmap_U,isRefAnat,isPopout);
end
%% threshold the bubble plot with chosen radius
thres_rad = 0.3;

radius = sqrt(stimcorr.^2 + motorcorr.^2);
[A,U_sorted] = sort(radius,'descend');
i_end = find(A>thres_rad,1,'last');
disp(i_end)
IX_passrad = U_sorted(1:i_end);

cIX_radthres = []; gIX_radthres = [];
for i = 1:length(IX_passrad),
    IX = find(gIX_in==IX_passrad(i));
    cIX_radthres = [cIX_radthres; cIX_in(IX)];
    gIX_radthres = [gIX_radthres; gIX_in(IX)];
end
%% Anat plot of only clusters > radius, same colormap
if isplotfig,
    isRefAnat = 1;
    isPopout = 1;
    figure
    DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_radthres,gIX_radthres,clrmap);
end

%% bubble plot in 2-D color, with radius drawn
if isplotfig,
    figure('Position',[500,500,300,250]);hold on;
    scatter(stimcorr,motorcorr,U_size,clrmap)
    xlabel('stimulus corr.');ylabel('motor corr.');
    theta = -1:0.01:pi/2;
    X = cos(theta)*thres_rad;
    Y = sin(theta)*thres_rad;
    plot(X,Y,'k--','Linewidth',1.5)
    axis equal
    % ylim([-0.15,0.4])
    % xlim([0,0.6])
    xlim([0,0.7]);
    ylim([-0.18,0.5]);
    set(gca,'YTick',[-0.2:0.2:0.4])
end

%%
thres_x = 0.4;% 0.1; for fish8 figure
[A,U_sorted] = sort(stimcorr,'ascend');
i_end = find(A<thres_x,1,'last');
IX_smallstim = U_sorted(1:i_end);
IX_plot = intersect(IX_smallstim,IX_passrad,'stable');
disp(length(IX_plot));

cIX_out = []; gIX_out = [];
for i = 1:length(IX_plot),
    IX = find(gIX_in==IX_plot(i));
    cIX_out = [cIX_out; cIX_in(IX)];
    gIX_out = [gIX_out; i*ones(size(IX))];
end
%%
cIX = cIX_out;
gIX = gIX_out;