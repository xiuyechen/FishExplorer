
[cIX_all,gIX_all,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[2,1]);


%%

i_lr = 2;
betas = Betas{i_lr,i_fish};
b1 = betas(:,1);
b2 = betas(:,2);
b3 = betas(:,3);


%%
caseflag = 2;
switch caseflag
    case 1
        X = b1;%b3;%b1;
        Y = b3;%b2;
        Z = b3;
        Xname = 'motor only (b1)';
        Yname = 'stim only (b3)';
    case 2
        X = b1;%b3;%b1;
        Y = sqrt(b2.^2+b3.^2);%b2;
        Z = b3;
        Xname = 'motor only (b1)';
        Yname = 'periodic';
end
% Xname = 'motor only (b1)';
% Yname = 'SMT (b2)';

A = Y;
numcell = length(X);
topN = round(0.05*numcell); % top 5% cutoff
[~,IX] = sort(A,'descend');
thresA = A(IX(topN));

% get min/max
x0 = min(X(IX(1:topN)));
x1 = max(X(IX(1:topN)));
y0 = min(Y(IX(1:topN)));
y1 = max(Y(IX(1:topN)));
z0 = min(Z(IX(1:topN)));
z1 = max(Z(IX(1:topN)));

IX_pass = find(A>=thresA);
IX_fail = find(A<thresA);

gIX_in = (1:length(X))';
%% make custom 2-D colormap
grid = Make4color2Dcolormap;
if false
    grid = MakeDiagonal2Dcolormap;
end

%% map data to 2D colormap, and cluster sizes
clrmap0 = MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
clrmap = clrmap0;
clrmap(IX_fail,:) = ones(length(IX_fail),3)*0.5;

if length(clrmap)>1000
    U_size = ones(size(X));
else % for actual clusters
    % get cluster size (number of cells in each cluster)
    gIX2 = SqueezeGroupIX(gIX_in);
    U = unique(gIX2);
    U_size = zeros(size(X));
    for i = 1:length(U)
        ix = find(gIX2 == U(i));
        U_size(i) = length(ix);
    end
end

%% bubble plot
h = figure('Position',[500,500,300,250]); hold on
scatter(X,Y,U_size,clrmap,'filled')
plot([x0,x1],[y0,y0],'k--');

xlabel(Xname);ylabel(Yname);
axis equal
xlim([-0.5,1]);
ylim([-0.5,1]);
%     set(gca,'YTick',-0.2:0.2:0.6);


%%
% scatter(b1(cIX),b2(cIX),10,'filled')
%%
c1 = SelectClusterRange(cIX,gIX,1:64);
[c2,g2] = SelectClusterRange(cIX,gIX,65:128);

scatter(X(c1),Y(c1),8,'r','filled');
scatter(X(c2),Y(c2),8,'g','filled');

%%
[cIX_hb,gIX_hb] = ScreenCellsWithMasks([219,220],cIX,gIX,MASKs,CellXYZ_norm,absIX);
[c2_hb,g2_hb] = SelectClusterRange(cIX_hb,gIX_hb,65:128);
scatter(X(c2_hb),Y(c2_hb),8,'k','filled');

%% find b1<0.2
thres = 0.15;
c_X = X(c);
ix_low = find(c_X<thres);
ix_high = find(c_X>thres);
g(ix_low) = 1;
g(ix_high) = 2;
% cIX_low = c2(ix_low);
% cIX_high = c2(ix_high);

setappdata(hfig,'clrmap_name','hsv_old');
% make figure
I = LoadCurrentFishForAnatPlot(hfig,c,g);
[h,im_full2] = DrawCellsOnAnat(I);

%%
c_Y = Y(c);
scatter(c_X(ix_low),c_Y(ix_low),4,'k','filled');
scatter(c_X(ix_high),c_Y(ix_high),4,'c','filled');

%% find x>0.4
ix_low = find(b1(c)<0.2);
ix_high = find(b1>0.5);
cL = c(ix_low);
cH = cIX_all(ix_high);

gL = ones(size(cL));
gH = 2*ones(size(cH));

cIX_plot = vertcat(cL,cH);
gIX_plot = vertcat(gL,gH);

setappdata(hfig,'clrmap_name','hsv_old');
% make figure
I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot);
[h,im_full2] = DrawCellsOnAnat(I);

%%
% cIX_plot = cIX;%c2_hb;
% gIX_plot = gIX;

% cont from C:\Users\xiuye\Dropbox\FishExplorer2\figure scripts\Fig5 MultiFish\Batch_PTintOMR_periodic_scatterplot_anatmap.m

        clr1 = [0.5,1,0.5];
        clr2 = [0.2,0.4,1];%[0.3,0.8,1];
%         clr_fail = [0.5,0.5,0.5];
        clr_int = [1,0.1,0];
        
%  IX_pass_2 = union(union(PassX_2{1},PassX_2{2}),union(PassY_2{1},PassY_2{2}));
 cIX_plot = [cIX_all(PassX_2{2});cIX_all(PassY_2{2});cIX_int];
 gIX_plot = [ones(size(cIX_all(PassX_2{2})));2*ones(size(cIX_all(PassY_2{2})));3*ones(size(cIX_int))];
%  clrmap = [repmat(clr1,length(gIX_int),1);...
%      repmat(clr2,length(gIX_int),1);...
%      repmat(clr_int,length(gIX_int),1)];
clrmap = [clr1;clr2;clr_int];
%%
figure('Position',[50,100,400,500]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',0);
setappdata(hfig,'stimrange',[1,2]);
UpdateIndices_Manual(hfig,cIX_plot,gIX_plot);
UpdateTimeIndex(hfig);

opts = [];
opts.clrmap = clrmap;
DrawTimeSeries(hfig,cIX_plot,gIX_plot,opts);
%%
cIX = cIX_plot;
gIX = gIX_plot;
%%
cIX_plot = cIX;
gIX_plot = gIX;
clrmap = [clr2;clr2];%clr_int;clr1];