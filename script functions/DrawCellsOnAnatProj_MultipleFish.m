function  [tot_image, dim_totimage] = DrawCellsOnAnatProj_MultipleFish(hfig,AllCentroids,range_fish,range_clus,isPopout)
anat_yx = getappdata(hfig,'anat_yx_norm');
anat_yz = getappdata(hfig,'anat_yz_norm');
anat_zx = getappdata(hfig,'anat_zx_norm');
k_zres = 2.5;
radius_xy = 3;
width_z = 5;
thickness_z = 3;

%%
CellXYZ = [];
gIX = [];

for i = 1:length(range_fish), 
    i_fish = range_fish(i);
    XYZ_thisfish = vertcat(AllCentroids{i_fish}.XYZn{range_clus{i}});
    CellXYZ = vertcat(CellXYZ,XYZ_thisfish);
    gIX = vertcat(gIX,i*ones(size(XYZ_thisfish,1),1));
end
cIX = (1:length(gIX))';
numK = length(range_fish);

%% load
clrmap = getappdata(hfig,'clrmap');
isShowMasks = getappdata(hfig,'isShowMasks');
if isShowMasks,
    MASKs = getappdata(hfig,'MASKs');
    Msk_IDs = getappdata(hfig,'Msk_IDs');   
    isShowMskOutline = getappdata(hfig,'isShowMskOutline');  
end
isWeighAlpha = getappdata(hfig,'isWeighAlpha');

%% Setup
[s1,s2] = size(cIX);
if s2>s1,
    cIX = cIX';
end
[s1,s2] = size(gIX);
if s2>s1,
    gIX = gIX';
end

% % down-sample
% if ~isPopout,
%     displaymax = 8000;
%     if length(cIX_abs) > displaymax,
%         skip = round(length(cIX_abs)/displaymax);
%         cIX_abs = cIX_abs(1:skip:end,:);
%         gIX = gIX(1:skip:end,:);
%     end
% end

% get numK
if exist('numK','var'),
    numK = double(max(numK,max(gIX)));
else
    numK = double(max(gIX));
end

% get colormap
if strcmp(clrmap,'jet'),
%     temp = zeros(numK,3);
%     temp(:,1) = linspace(1,0,numK);
%     temp(:,2) = linspace(0,0,numK);
%     temp(:,3) = linspace(0,0,numK);
%     cmap = temp;
%     cmap = flipud(gray(double(numK)));
    cmap = flipud(jet(numK));
else % 'hsv'
    n = round(numK*1.1);
    cmap = hsv(max(1,n));
end

% set up anatomical images
anat_YX = anat_yx/4;
anat_YZ = anat_yz/4;
anat_ZX = anat_zx/4;
anat_ZX = flipud(anat_ZX);
dimv_yx = size(anat_YX);
dimv_yz = size(anat_YZ);
dimv_zx = size(anat_ZX);

% adjust size for YZ and ZX view, and combine all 3 view into total-image
% (including white borders as indicated in dim_totimage)
anat_yz2=zeros(dimv_yz(1),dimv_yz(2)*k_zres,3);
anat_zx2=zeros(dimv_zx(1)*k_zres,dimv_zx(2),3);
dim_totimage = [dimv_yx(1)+dimv_zx(1)*k_zres+10,dimv_yx(2)+dimv_yz(2)*k_zres+10,3];
tot_image=ones(dim_totimage);

% create circlular mask (~linearized indices), to draw cells on image (later)
circle=makeDisk2(radius_xy,radius_xy*2+1); % make mask of filled circle % (7,15)
mask = zeros(dimv_yx(1),dimv_yx(2));
mask(1:radius_xy*2+1,1:radius_xy*2+1) = circle;
ix = find(mask);
cix = sub2ind([dimv_yx(1),dimv_yx(2)],radius_xy+1,radius_xy+1);% 8
circle_inds = ix - cix;

% set transparancy
alpha_max = 0.3-min(length(cIX)/1000/100,0.1);
if isWeighAlpha,    
    wIX = getappdata(hfig,'wIX');
    alpha = wIX*alpha_max;
else
    alpha = ones(size(cIX))*alpha_max;
end

%% Draw each cell on map (all 3 views)
% set up specialized indices
square=ones(width_z,thickness_z); % make mask of filled circle % (7,15)
mask = zeros(dimv_yz(1),dimv_yz(2));
mask(1:width_z,1:thickness_z) = square;
ix = find(mask);
cix = sub2ind([dimv_yz(1),dimv_yz(2)],round((width_z+1)/2),round((thickness_z+1)/2));
yzplane_inds = ix - cix;
% yzplane_inds = -radius_z:radius_z;

square=ones(thickness_z,width_z); % make mask of filled circle % (7,15)
mask = zeros(dimv_zx(1),dimv_zx(2));
mask(1:thickness_z,1:width_z) = square;
ix = find(mask);
cix = sub2ind([dimv_zx(1),dimv_zx(2)],round((thickness_z+1)/2),round((width_z+1)/2));
zxplane_inds = ix - cix;
% zxplane_inds = -width_z*dimv_zx(1):dimv_zx(1):width_z*dimv_zx(1);

% main loop
% cIX_abs = absIX(cIX_abs);
for j=1:length(cIX)        
    if ~isempty(CellXYZ(cIX(j),1)),
        ix = gIX(j);
                
        %% Y-X        
        cinds=(CellXYZ(cIX(j),2)-1)*dimv_yx(1)+CellXYZ(cIX(j),1); % linear pixel index, faster equivalent of:
        %     cinds = sub2ind([dim_y,dim_x],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).center(2));
        labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yx(1)*dimv_yx(2)); % within bounds
        ixs = cinds+circle_inds(labelinds);
        anat_YX(ixs) = min(1,(anat_YX(ixs)+cmap(ix,1)*alpha(j)))*0.8 + anat_YX(ixs)*0.2; % R
%         anat_YX(ixs) = cmap(ix,1)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % R
        ixs = cinds+circle_inds(labelinds) + dimv_yx(1)*dimv_yx(2);
        anat_YX(ixs) = min(1,(anat_YX(ixs)+cmap(ix,2)*alpha(j)))*0.8 + anat_YX(ixs)*0.2; % G
%         anat_YX(ixs) = cmap(ix,2)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % G
        ixs = cinds+circle_inds(labelinds) + dimv_yx(1)*dimv_yx(2)*2;
        anat_YX(ixs) = min(1,(anat_YX(ixs)+cmap(ix,3)*alpha(j)))*0.8 + anat_YX(ixs)*0.2; % B
%         anat_YX(ixs) = cmap(ix,3)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % B

        %% Y-Z
        z_alpha(j) = alpha(j);%/2;        
        cinds=(CellXYZ(cIX(j),3)-1)*dimv_yz(1)+CellXYZ(cIX(j),1); % linear pixel index, faster equivalent of:
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        labelinds=find((cinds+yzplane_inds)>0 & (cinds+yzplane_inds)<=dimv_yz(1)*dimv_yz(2));
        ixs = cinds+yzplane_inds(labelinds);
        anat_YZ(ixs) = min(1,(anat_YZ(ixs)+cmap(ix,1)*z_alpha(j)))*0.8 + anat_YZ(ixs)*0.2; % R
%         anat_YZ(ixs) = cmap(ix,1)*z_alpha(j) + anat_YZ(ixs)*(1-z_alpha(j)); % R
        ixs = cinds+yzplane_inds(labelinds) + dimv_yz(1)*dimv_yz(2);
        anat_YZ(ixs) = min(1,(anat_YZ(ixs)+cmap(ix,2)*z_alpha(j)))*0.8 + anat_YZ(ixs)*0.2; % G
%         anat_YZ(ixs) = cmap(ix,2)*z_alpha(j) + anat_YZ(ixs)*(1-z_alpha(j)); % G
        ixs = cinds+yzplane_inds(labelinds) + dimv_yz(1)*dimv_yz(2)*2;
        anat_YZ(ixs) = min(1,(anat_YZ(ixs)+cmap(ix,3)*z_alpha(j)))*0.8 + anat_YZ(ixs)*0.2; % B
%         anat_YZ(ixs) = cmap(ix,3)*z_alpha(j) + anat_YZ(ixs)*(1-z_alpha(j)); % B
                
        %% Z-X
        z_alpha(j) = alpha(j);%/2;        
        cinds=(CellXYZ(cIX(j),2)-1)*dimv_zx(1) +(CellXYZ(cIX(j),3)); % linear pixel index, faster equivalent of:
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        labelinds=find((cinds+zxplane_inds)>0 & (cinds+zxplane_inds)<=dimv_zx(1)*dimv_zx(2));
        ixs = cinds+zxplane_inds(labelinds);
        anat_ZX(ixs) = min(1,(anat_ZX(ixs)+cmap(ix,1)*z_alpha(j)))*0.8 + anat_ZX(ixs)*0.2; % R
%         anat_ZX(ixs) = cmap(ix,1)*z_alpha(j) + anat_ZX(ixs)*(1-z_alpha(j)); % R
        ixs = cinds+zxplane_inds(labelinds) + dimv_zx(1)*dimv_zx(2);
        anat_ZX(ixs) = min(1,(anat_ZX(ixs)+cmap(ix,2)*z_alpha(j)))*0.8 + anat_ZX(ixs)*0.2; % G
%         anat_ZX(ixs) = cmap(ix,2)*z_alpha(j) + anat_ZX(ixs)*(1-z_alpha(j)); % G
        ixs = cinds+zxplane_inds(labelinds) + dimv_zx(1)*dimv_zx(2)*2;
        anat_ZX(ixs) = min(1,(anat_ZX(ixs)+cmap(ix,3)*z_alpha(j)))*0.8 + anat_ZX(ixs)*0.2; % B
%         anat_ZX(ixs) = cmap(ix,3)*z_alpha(j) + anat_ZX(ixs)*(1-z_alpha(j)); % B
        
    end
end

% draw masks on anat projections
if isShowMasks,  
    nMasks = length(Msk_IDs);
    % set params
    cmap2 = jet(nMasks);
    msk_alpha(j) = 0.25;
    white_alpha(j) = 0.05;
    outline_radius = 5;
    
    for i = 1:nMasks,
        msk = full(MASKs.MaskDatabase(:,Msk_IDs(i)));
        mask_3D = reshape(msk, [MASKs.height, MASKs.width, MASKs.Zs]);
        
        % Y-X
        masks_XY = max(mask_3D,[],3);
        masks_fit = logical(imresize(masks_XY,[size(anat_yx,1),size(anat_yx,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), strel('disk',outline_radius));
        end
        ixs = find(masks_fit);
        anat_YX(ixs) = (cmap2(i,1)*msk_alpha(j) + anat_YX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % R;
        ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2);
        anat_YX(ixs) = (cmap2(i,2)*msk_alpha(j) + anat_YX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % G;
        ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2)*2;
        anat_YX(ixs) = (cmap2(i,3)*msk_alpha(j) + anat_YX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % B;
        
        % Y-Z
        masks_YZ = squeeze(max(mask_3D,[],2));
        masks_fit = logical(imresize(masks_YZ,[size(anat_yz,1),size(anat_yz,2)]));
        if isShowMskOutline,
%             masks_fit = edge(masks_fit);
            masks_fit = imdilate(edge(masks_fit), strel('line',outline_radius,90));
        end
        ixs = find(masks_fit);
        anat_YZ(ixs) = (cmap2(i,1)*msk_alpha(j) + anat_YZ(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % R;
        ixs = find(masks_fit) + dimv_yz(1)*dimv_yz(2);
        anat_YZ(ixs) = (cmap2(i,2)*msk_alpha(j) + anat_YZ(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % G;
        ixs = find(masks_fit) + dimv_yz(1)*dimv_yz(2)*2;
        anat_YZ(ixs) = (cmap2(i,3)*msk_alpha(j) + anat_YZ(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % B;
        
        % Z-X
        masks_ZX = squeeze(max(mask_3D,[],1))';  % notice the transpose
        masks_fit = logical(imresize(masks_ZX,[size(anat_zx,1),size(anat_zx,2)]));
        if isShowMskOutline,
%             masks_fit = edge(masks_fit);
            masks_fit = imdilate(edge(masks_fit), strel('line',outline_radius,0));
        end
        ixs = find(masks_fit);
        anat_ZX(ixs) = (cmap2(i,1)*msk_alpha(j) + anat_ZX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % R;
        ixs = find(masks_fit) + dimv_zx(1)*dimv_zx(2);
        anat_ZX(ixs) = (cmap2(i,2)*msk_alpha(j) + anat_ZX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % G;
        ixs = find(masks_fit) + dimv_zx(1)*dimv_zx(2)*2;
        anat_ZX(ixs) = (cmap2(i,3)*msk_alpha(j) + anat_ZX(masks_fit)*(1-msk_alpha(j)))*(1-white_alpha(j)) + white_alpha(j); % B;
    end
end

% rescale (low-res) z dimension
for k=1:3
    anat_yz2(:,:,1) = imresize(anat_YZ(:,:,1), [dimv_yz(1), dimv_yz(2)*k_zres],'nearest');
    anat_yz2(:,:,2) = imresize(anat_YZ(:,:,2), [dimv_yz(1), dimv_yz(2)*k_zres],'nearest');
    anat_yz2(:,:,3) = imresize(anat_YZ(:,:,3), [dimv_yz(1), dimv_yz(2)*k_zres],'nearest');
    
    anat_zx2(:,:,1) = imresize(anat_ZX(:,:,1), [dimv_zx(1)*k_zres, dimv_zx(2)],'nearest');
    anat_zx2(:,:,2) = imresize(anat_ZX(:,:,2), [dimv_zx(1)*k_zres, dimv_zx(2)],'nearest');
    anat_zx2(:,:,3) = imresize(anat_ZX(:,:,3), [dimv_zx(1)*k_zres, dimv_zx(2)],'nearest');
end

tot_image(dimv_zx(1)*k_zres+11:end,1:dimv_yx(2),:) = anat_YX;
tot_image(dimv_zx(1)*k_zres+11:end,dimv_yx(2)+11:end,:) = anat_yz2;
tot_image(1:dimv_zx(1)*k_zres,1:dimv_zx(2),:) = flipud(anat_zx2);

tot_image(tot_image(:)>1) = 1;
tot_image(tot_image(:)<0) = 0;

image(tot_image);
axis image;axis off

end

function out = makeDisk2(radius, dim)
center=floor(dim/2)+1;
out=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<=radius
            out(x,y)=1;
        end
    end
end

end