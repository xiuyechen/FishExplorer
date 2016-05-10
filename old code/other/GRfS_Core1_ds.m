%%% Grow ROI's from Seeds
if projnum == 1, % for HBO
    thres_minarea = 10; 
    thres_area = 30;
else % if projnum ==2, % for spatial PT
    thres_minarea = 7;
    thres_area = 20;
end
thres_seed = 0.35;
thres_nb = 0.3;
% mfsize = 1; % for spatial 2D median filter on raw correlation stacks

%% generate dF/F over all pixels
disp('normalize IM...')
Y = double(reshape(IM,[s1*s2 nFrames])); % time-course
Y_mean = nanmean(Y,2)*ones(1,nFrames); % should be no NaN's
dFF_IM = (Y-Y_mean)./Y_mean;
clear Y, clear Y_mean;

%% calculate correlation with neighboring pixels % takes ~15sec for 160x100 px
IM1 = reshape(dFF_IM,[s1,s2,nFrames]);

% corr with 3x3 neighboring pixels calculated, and averaged.
translation = [-1,-1;-1,0;-1,1;0,-1;0,1;1,-1;1,0;1,1]; % 8 directions
% [I,J]=ind2sub([3,3],1:9); ... handcoding is easier for 3x3
temp = zeros(s1,s2,8);
for i=1:8,
    IM2 = imtranslate(IM1,translation(i,:));
    A = reshape(IM1,[s1*s2, nFrames])'; 
    B = reshape(IM2,[s1*s2, nFrames])'; 
    
    %% below is the faster ~equivalent (super close) of diag(corr(A,B))
    An=bsxfun(@minus,A,mean(IM1(:),1));
    Bn=bsxfun(@minus,B,mean(IM2(:),1));
    An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
    Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
    C=sum(An.*Bn,1);
    disp(['translation #' num2str(i)]);
    temp(:,:,i) = reshape(C,[s1,s2]);
end
im_corr_nb(:,:,i_plane) = mean(temp,3);
clear IM2;

%% ITERATION:
% initialize
disp('find ROIs')
ROI = [];
ROI(300).pxlist = []; % 300 not capping
ROI(300).area = [];

BW = zeros(s1,s2); % mask to store new ROI
ROIcount = 0;

% find 1st seed
m = im_corr_nb(:,:,i_plane);
[C,IX] = max(m(:));
[I,J] = ind2sub([s1,s2],IX);

while m(I,J)>thres_seed, % each seed
    if (I==1 || I==s1 || J==1 || J==s2),% edge of image
        % find next seed and skip to next loop
        m(I,J) = NaN;
        [C,IX] = max(m(:));
        [I,J] = ind2sub([s1,s2],IX);
        continue;
    end
    % add seed position to mask
    BW(I,J) = 1;

    % grow ROI: expand mask and examine next neighbors
    IX2 = 0; area = 0; % initialize to dummy
    while ~isempty(IX2) && area<=(thres_area-4), % 4 to grow this round
        BW0 = imdilate(BW,[0 1 0; 1 1 1; 0 1 0]);% 4 neighbors
        BW0(find(BW))=0;
        IX1 = find(BW0);
        IX2 = find(m(IX1)>thres_nb);
        if ~isempty(IX2),
            BW(IX1(IX2)) = 1;
            area = length(find(BW));
        end
    end % finished expanding this seed

    % save this ROI if bigger than thres
    area = length(find(BW));
    if area >= thres_minarea,% = 3 for ds 0.2
        ROIcount = ROIcount+1;
        ROI(ROIcount).pxlist = find(BW);
        ROI(ROIcount).area = area;
    end
    % reset
    m(find(BW)) = NaN;
    BW = zeros(s1,s2);
    % look for next seed
    [C,IX] = max(m(:));
    [I,J] = ind2sub([s1,s2],IX);
end

% no seeds above thres_nb found. 
ROI(ROIcount+1:end) = [];

%% save fluo-traces
for i = 1:ROIcount,
    A = reshape(IM, s1*s2, nFrames); % rows: pixels within ROI; cols: fluo-trace
    B = A(ROI(i).pxlist,:);
    fluotrace1D = sum(B,1);
    ROI(i).fluo_trace = fluotrace1D / ROI(i).area;% zscore normalize to mean = 0, SD = 1
    ROI(i).fluoZ = zscore(fluotrace1D);% zscore normalize to mean = 0, SD = 1
end

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
