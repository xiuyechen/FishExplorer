

%%

cIX_in = getappdata(hfig,'cIX');
gIX_in = getappdata(hfig,'gIX_betas');

%% get coefficients from multiple linear regression
betas = getappdata(hfig,'betas');
stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
% figure;scatter(stimcorr,motorcorr)%scatter(motorcorr,stimcorr)

%% draw custom 2-D colormap illustration (square)
res = 100;

huex = 120/360; % 120/360 for green/magenta
%huey = 280/360;
satmin = 0;
pw = 1; % 
grid = zeros(res,res,3);

% grid(:,:,3) = repmat(grad,res,1)';
% grid(:,:,1) = 0.5*repmat(rev_grad',1,res)'+0.5*repmat(rev_grad,res,1)';
% grid(:,:,2) = repmat(grad',1,res)';

for i = 1:res
    for j = 1:res
%         grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(i*colorx/res)';
%         grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(j*colory/res)';
        % hue(i,j) = (i*huex + j*huey)/(2*res); % average the hue
        val(i,j) = min((sqrt(i^2+j^2)/res)^pw,1);
            if i>j
                hue(i,j) = huex;
                sat(i,j) = max(satmin,(i - j)/res);
            else
                hue(i,j) = huex+0.5;
                sat(i,j) = max(satmin,(j - i)/res);
            end
        grid(i,j,:) = hsv2rgb([hue(i,j), sat(i,j), val(i,j)]);
        % grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(j*colory/res)';

    end
end

clrmap_2D = reshape(grid,res*res,3);

figure;imagesc(grid)
axis xy
axis off
axis equal

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
figure('Position',[500,500,300,250]);
scatter(stimcorr,motorcorr,U_size,clrmap)
xlabel('stimulus corr.');ylabel('motor corr.');
axis equal
xlim([0,0.7]);
ylim([-0.22,0.6]);
% set(gca,'YTick',-0.2:0.2:0.6);

%% Anat plot with custom colormap
isRefAnat = 1;
isPopout = 1;
figure
DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_in,gIX_in,clrmap);
% DrawCellsOnAnatProj_othercolor(hfig,cIX_in,gIX_in,cmap_U,isRefAnat,isPopout);

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
isRefAnat = 1;
isPopout = 1;
figure
DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_radthres,gIX_radthres,clrmap);

%% bubble plot in 2-D color, with radius drawn
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
ylim([-0.22,0.6]);
set(gca,'YTick',-0.2:0.2:0.6);

%%
thres_x = 0.3;% 0.1; for fish8 figure
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

%% draw custom 3-D colormap illustration (square)
res = 20;
grad = linspace(0,1,res);
rev_grad = linspace(1,0,res);

colorx = [1,1,0];
colory = [0,0,1];


hue1 = 60/360;
hue2 = 180/360;
hue3 = 300/360;
%huey = 280/360;
satmin = 0.3;
pw = 5;
grid = zeros(res,res,res,3);

% grid(:,:,3) = repmat(grad,res,1)';
% grid(:,:,1) = 0.5*repmat(rev_grad',1,res)'+0.5*repmat(rev_grad,res,1)';
% grid(:,:,2) = repmat(grad',1,res)';


for i = 1:res
%      tic
    for j = 1:res
        for k = 1:res
           
%         grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(i*colorx/res)';
%         grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(j*colory/res)';
        % hue(i,j) = (i*huex + j*huey)/(2*res); % average the hue
        val(i,j,k) = min((sqrt(i^2+j^2+k^2)/res)^pw,1);
        [M,I] = max([i,j,k]);
        switch I
            case 1
                hue(i,j,k,:) = hue1;
                sat(i,j,k,:) = max(satmin,(i - (j + k)/2)/res);
            case 2
                hue(i,j,k,:) = hue2;
                sat(i,j,k,:) = max(satmin,(j - (i + k)/2)/res);
            case 3
                hue(i,j,k,:) = hue3;
                sat(i,j,k,:) = max(satmin,(k - (i + i)/2)/res);
        end
        
        color = hsv2rgb([hue(i,j,k), sat(i,j,k), val(i,j,k)]);
        grid(i,j,k,:) = color;

        % grid(i,j,:) = squeeze(grid(i,j,:))+squeeze(j*colory/res)';
        
        end
    end
%     toc
end
%%
clrmap_3D = reshape(grid,res*res*res,3);

n = 10000;
xidx = randi(res,n,1);
yidx = randi(res,n,1);
zidx = randi(res,n,1);

figure; hold on;
for i = 1:n
    color(i,:) = squeeze(grid(xidx(i),yidx(i),zidx(i),:))';
end
scatter3(xidx,yidx,zidx,20,color,'filled');
view(120,30);
box on

%% get new gIX with matching custom colormap 'cmap_U'
gIX_x = round((stimcorr-min(stimcorr))/(max(stimcorr)-min(stimcorr))*(res-1))+1;
gIX_y = round((motorcorr-min(motorcorr))/(max(motorcorr)-min(motorcorr))*(res-1))+1;
gIX_z = repmat(14,length(gIX_x),1);


gIX_old = gIX_in;
U = unique(gIX_old);
U_size = zeros(size(gIX_x));
clrmap = zeros(length(gIX_x),3);
for i = 1:length(U);
    ix = find(gIX_old == U(i));
    U_size(i) = length(ix);
    ix = sub2ind([res,res,res],gIX_z(U(i))',gIX_y(U(i))',gIX_x(U(i))');
    clrmap(i,:) = clrmap_3D(ix,:);
end

%% Anat plot with custom colormap
isRefAnat = 1;
isPopout = 1;
figure
DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_in,gIX_in,clrmap);
% DrawCellsOnAnatProj_othercolor(hfig,cIX_in,gIX_in,cmap_U,isRefAnat,isPopout);