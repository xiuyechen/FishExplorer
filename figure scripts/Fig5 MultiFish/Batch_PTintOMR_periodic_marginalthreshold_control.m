% This script is for the control for the scatterplot in the last fig, taking a 
% 'slice' of leniently defined PT&OMR intersection cells are plotting them
% in the scatter plot, with anatomical location (anterior Hb vs posterior)
% labeled. Continues from middle (~120 lines) of
% 'Batch_PTintOMR_periodic_thresholdsweep_sandbox.m'.
% Not yet polished for presentation.

%%
i_prct_count = 9;
cIX1 = M_cIX{i_prct_count,1};
cIX2 = M_cIX{i_prct_count,2};

[cIX_int_low,ix] = intersect(cIX1,cIX2);
%         gIX_int_low = ones(size(cIX_int));% not used here but for other plotting

i_prct_count = 10;
cIX1 = M_cIX{i_prct_count,1};
cIX2 = M_cIX{i_prct_count,2};

[cIX_int_high,ix] = intersect(cIX1,cIX2);
%         gIX_int = 2*ones(size(cIX_int));% not used here but for other plotting

cIX_int = setdiff(cIX_int_high,cIX_int_low);
gIX_int = ones(size(cIX_int));

%% draw anat map
setappdata(hfig,'clrmap_name','jet');
I = LoadCurrentFishForAnatPlot(hfig,cIX_int,gIX_int);%,clrmap);
[h,~,im] = DrawCellsOnAnat(I);

%% below: make bubble-plot

load(fullfile(outputDir,'4D_SM_stimrangePTOMR_betas.mat'));
betas = Betas{i_lr,i_fish};
b1 = betas(:,1);
b2 = betas(:,2);
b3 = betas(:,3);
%%
% set up plot dimensions
caseflag = 2;
switch caseflag
    case 1
        X = b1;
        Y = b2;
        Xname = 'motor only (b1)';
        Yname = 'SMT (b2)';
    case 2
        X = b1;%b3;%b1;
        Y = sqrt(b2.^2+b3.^2);%b2;
        Xname = 'motor only (b1)';
        Yname = 'periodic';
end

numcell = length(b1);

A = X;
topN = length(M_cIX{2});%length(cIX_int);%round(0.01*numcell); % top 5% cutoff
[~,IX] = sort(A,'descend');
thresA = A(IX(topN));

B = Y;
topN = length(M_cIX{2});%length(cIX_int);%round(0.01*numcell); % top 5% cutoff
[~,IX] = sort(B,'descend');
thresB = min(Y(cIX_int));%B(IX(topN));

IX_passX = setdiff(find(A>=thresA),find(B>=thresB));
IX_passY = find(B>=thresB);%setdiff(find(B>=thresB),find(A>=thresA));
IX_pass = union(find(A>=thresA),find(B>=thresB));
IX_fail = intersect(find(A<thresA),find(B<thresB));%find(A<thresA);

% get min/max
x0 = min(X(IX_pass));
x1 = max(X(IX_pass));
y0 = min(Y(IX_pass));
y1 = max(Y(IX_pass));

gIX_in = (1:length(X))';

PassX_2{i_lr} = IX_passX;
PassY_2{i_lr} = IX_passY;
%% make colormap

clr1 = [0.3,0.8,0.2];
clr2 = [0.1,0.3,1];%[0.3,0.8,1];
clr_fail = [0.5,0.5,0.5];
clr_int = [1,0.1,0];
clrmap = ones(numcell,3);%MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
clrmap(IX_passX,:) = clr1.*ones(length(IX_passX),3);
clrmap(IX_passY,:) = clr2.*ones(length(IX_passY),3);
clrmap(IX_fail,:) = clr_fail.*ones(length(IX_fail),3);

%% bubble plot
h = figure('Position',[500,100,300,250]); hold on
U_size = ones(size(X));
scatter(X,Y,U_size,clrmap,'filled')
plot([x0,x1],[thresB,thresB],'k--');
plot([thresA,thresA],[y0,y1],'k--');
%         plot([x0,x1],[y0,y0],'k--');

scatter(X(IX_passX),Y(IX_passX),2,clr1);%,'filled');
scatter(X(IX_passY),Y(IX_passY),2,clr2);%,'filled');

xlabel(Xname);ylabel(Yname);
axis equal

xlim([-0.4,0.8]);
ylim([0,1]);

scatter(X(cIX_int),Y(cIX_int),2,clr_int);%[1,0.5,0.5]);%,'filled');

%% highlight anterior hindbrain (Rh1&2) cells (Rh1 219; Rh2 220;)
MASKs = getappdata(hfig,'MASKs');
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
absIX = getappdata(hfig,'absIX');
[c_ahb,g_ahb] = ScreenCellsWithMasks([219,220],cIX_int,gIX_int,MASKs,CellXYZ_norm,absIX);
clr_ahb = [1,0.8,0];
scatter(X(c_ahb),Y(c_ahb),5,clr_ahb);

[c_phb,g_phb] = ScreenCellsWithMasks([221,222,223,224],cIX_int,gIX_int,MASKs,CellXYZ_norm,absIX);
clr_phb = [0,1,1];
scatter(X(c_phb),Y(c_phb),5,clr_phb);


