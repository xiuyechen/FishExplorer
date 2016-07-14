% Set params
thres_minsize = 10;

thres_seed = 0.25;
thres_grow = 0.25;
thres_cap = 0.3;

%% generate dF/F over all pixels
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');

gIX = f.HierClus(M,gIX);
C = FindCentroid(hfig);
nPixels = size(C,1);

%% Calculate correlation distance between all clusters
temp = pdist(C,'correlation');
Dist = squareform(temp);
for i = 1:nPixels,
    Dist(i,i) = NaN;
end

%% ITERATION:
% initialize
disp('find ROIs')
ROI = [];
ROI(nPixels).pxlist = []; % not capping
ROI(nPixels).numcell = [];

BW = zeros(nPixels,1); % analogy ~ mask to store new ROI
ROIcount = 0;

% find 1st seed
[~,ix] = min(Dist(:));
[I,J] = ind2sub(size(Dist),ix);

%%
while Dist(I,J)<thres_seed, % each seed
    %% add seed position to mask
    BW(I) = 1;
    BW(J) = 1;

    % find new core of seed based on cells in clusters
    IX1 = find(gIX==I);
    IX2 = find(gIX==J);
    M_core = vertcat(M(IX1,:),M(IX2,:));
    [~,newCore] = kmeans(M_core,1,'distance','correlation');

    % find next closest cluster
    D = 1-corr(newCore',C'); % distance matrix between new core and the original cluster means
    IX = find(BW==0);
    [a,ix] = min(D(IX));
    
    %% grow ROI: expand mask and examine next neighbors
    while true
        % find new core based on cells in clusters
        list = find(BW);
        IX = [];
        for i_gIX = 1:length(list),
            IX = [IX;find(gIX==list(i_gIX))];
        end
        M_core = M(IX,:);
        [~,newCore] = kmeans(M_core,1,'distance','correlation');
        
        % find next closest cluster
        D = 1-corr(newCore',C'); % distance matrix between new core and the original cluster means
        IX1 = find(BW==0);
        [a,IX2] = min(D(IX1));
        
        BW_pretend = BW;
        BW_pretend(IX1(IX2)) = 1;
        IX = find(BW_pretend);
        maxdist = max(pdist(C(IX,:),'correlation'));    
        
        if a<thres_grow && maxdist<thres_cap,
            BW(IX1(IX2)) = 1;
        else
            break; % finished expanding this seed
        end
    end        
    
    % save this ROI if bigger than thres    
    list = find(BW);
    IX = [];
    for i_gIX = 1:length(list),
        IX = [IX;find(gIX==list(i_gIX))];
    end
    numcell = length(IX);
    if numcell >= thres_minsize,
        ROIcount = ROIcount+1;
        ROI(ROIcount).pxlist = find(BW);
        ROI(ROIcount).numcell = numcell;
    end
    
    % reset
    Dist(list,:) = NaN;
    Dist(:,list) = NaN;
    C(list,:) = NaN;
    BW = zeros(nPixels,1); 
    
    % look for next seed
    [~,ix] = min(Dist(:));
    [I,J] = ind2sub(size(Dist),ix);
end
 
ROI(ROIcount+1:end) = [];

%% Merge accordingly
for i = 1:ROIcount,
   pxlist = ROI(i).pxlist;
   for j = 2:length(pxlist),
      gIX(gIX==pxlist(j)) = pxlist(1);
   end
end
gIX = SqueezeGroupIX(gIX);

f.UpdateIndices(hfig,cIX,gIX);
SaveCluster_Direct(hfig,cIX,gIX,'MergefromSeed0.25',1);
UpdateClustersGUI(hfig);
f.RefreshFigure(hfig);

%% Regression with the centroid of each cluster, round 2
disp('auto-reg-clus, round 2');
[cIX,gIX] = AllCentroidRegression_direct(hfig);
f.UpdateIndices(hfig,cIX,gIX);
f.RefreshFigure(hfig);

%%
thres_size = getappdata(hfig,'thres_size');
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_size,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX, numU] = SqueezeGroupIX(gIX);
f.UpdateIndices(hfig,cIX,gIX);

gIX = f.HierClus(M,gIX);
f.UpdateIndices(hfig,cIX,gIX);
SaveCluster_Direct(hfig,cIX,gIX,'fromSeed_thres',3);

f.RefreshFigure(hfig);

%% Visualize all ROI's in ds resolution (optional background im_corr_nb)
% cmap = jet(ROIcount);
if 0
    disp('Draw ROIs alone')
    cmap = rand(ROIcount,3);
    clrim = zeros(s1,s2,3);
    for i = 1:ROIcount,
        im = zeros(s1,s2);
        im(ROI(i).pxlist) = 1;
        for k = 1:3,
            clrim(:,:,k) = clrim(:,:,k)+im * cmap(i,k);
        end
    end
    % add background
    % for k = 1:3,
    %     clrim(:,:,k) = im_corr_nb(:,:,i_plane)*0.5 + clrim(:,:,k)*0.5;
    % end
    figure; imagesc(clrim)
    set(gca,'dataAspectRatio',[1 1 1],'YDir','normal'); axis off
end
%% Draw on top of anatomy, full res, smudged
disp('Draw on anatomy')
im_anat = mat2gray(anatomy_stack(:,:,i_plane));
height = size(anatomy_stack,1);
width = size(anatomy_stack,2);

cmap = rand(ROIcount,3);
% cmap = jet(ROIcount);
clrim = zeros(height,width,3);
GaussFilter = fspecial('gaussian',[8,8],8);

figure; hold on
for k = 1:3,
        clrim(:,:,k) = im_anat;
end
for i = 1:ROIcount,
	BW = zeros(s1,s2);
    BW(ROI(i).pxlist) = 1;
    BW2 = imresize(BW,[height,width]);
%     BW2 = imfilter(BW2,GaussFilter);
    BW2(BW2<0)=0;
%     BW2 = mat2gray(BW2);
    for k = 1:3,
        clrim(:,:,k) = clrim(:,:,k) + BW2*cmap(i,k);
    end
end
clrim(clrim<0)=0;clrim(clrim>1)=1;

% draw
imagesc(clrim);
set(gca,'dataAspectRatio',[1 1 1],'YDir','normal'); axis off

%% pool (this plane)
im_ROI_stack(:,:,:,i_plane) = clrim;
ROI_stack{i_plane} = ROI;

% NOT CUTTING UP ROI'S LIKE RUBEN'S PAPER
% STATS = regionprops(BW, 'Orientation','Centroid','MajorAxisLength');
