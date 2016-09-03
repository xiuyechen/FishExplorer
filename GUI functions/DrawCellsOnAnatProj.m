function  [tot_image, dim_totimage] = DrawCellsOnAnatProj(hfig,isRefAnat,isPopout,cIX_plot,gIX_plot,clrmap)
%% load
if ~isRefAnat, % raw images
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres_ratio = 20; % 8/0.406 = 19.704
    radius_xy = 5;%7;
    width_z = 10;
    thickness_z = 1;
else % registered to ZBrain
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_yx = getappdata(hfig,'anat_yx_norm');
    anat_yz = getappdata(hfig,'anat_yz_norm');
    anat_zx = getappdata(hfig,'anat_zx_norm');
    k_zres_ratio = 2.5; % 2/0.798 = 2.506
    radius_xy = 3;
    width_z = 5;
    thickness_z = 3;    
    z_range_ventral = round(min(CellXYZ(:,3))*k_zres_ratio);
        
    % fish outline
    isShowFishOutline = getappdata(hfig,'isShowFishOutline');
    if isShowFishOutline,
        FishOutline = getappdata(hfig,'FishOutline');
        % option with some background showing
%         anat_yx = 0.5*anat_yx + 0.5*repmat(FishOutline.outline_XY,[1,1,3]);
%         anat_yz = 0.5*anat_yz + 0.5*repmat(FishOutline.outline_YZ,[1,1,3]);
%         anat_zx = 0.5*anat_zx + 0.5*repmat(FishOutline.outline_ZX,[1,1,3]);
        
        anat_yx = repmat(FishOutline.outline_XY,[1,1,3]);
        anat_yz = repmat(FishOutline.outline_YZ,[1,1,3]);
        anat_zx = repmat(FishOutline.outline_ZX,[1,1,3]);
    end
end

absIX = getappdata(hfig,'absIX');
if exist('cIX_plot','var'),
    cIX = cIX_plot;
else
    cIX = getappdata(hfig,'cIX');
end
if exist('gIX_plot','var'),
    gIX = gIX_plot;
else
    gIX = getappdata(hfig,'gIX');
end

if ~exist('clrmap','var'),
    % get colormap
    clrmap_name = getappdata(hfig,'clrmap_name');
    numK = getappdata(hfig,'numK');
    numK = double(max(numK,max(gIX)));% sadly this is not always true
    clrmap = GetColormap(clrmap_name,numK);
end

isShowMasks = getappdata(hfig,'isShowMasks');
if isShowMasks,
    isShowMskOutline = getappdata(hfig,'isShowMskOutline');  
    MASKs = getappdata(hfig,'MASKs');
    Msk_IDs = getappdata(hfig,'Msk_IDs');   
    if Msk_IDs == 0,
        newMask = getappdata(hfig,'newMask');
    end    
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

% down-sample
if ~isPopout,
    displaymax = 8000;
    if length(cIX) > displaymax,
        skip = round(length(cIX)/displaymax);
        cIX = cIX(1:skip:end,:);
        gIX = gIX(1:skip:end,:);
    end
end

% set up anatomical images
% anat_YX = zeros(size(anat_yx));
% anat_YZ = zeros(size(anat_yz));
% anat_ZX = zeros(size(anat_zx));
anat_YX = anat_yx/4;
anat_YZ = anat_yz/4;
anat_ZX = anat_zx/4;
anat_ZX = flipud(anat_ZX);
dimv_yx = size(anat_YX);
dimv_yz = size(anat_YZ);
dimv_zx = size(anat_ZX);

% adjust size for YZ and ZX view, and combine all 3 view into total-image
% (including white borders as indicated in dim_totimage)

anat_yz2=zeros(dimv_yz(1),dimv_yz(2)*k_zres_ratio,3);
anat_zx2=zeros(dimv_zx(1)*k_zres_ratio,dimv_zx(2),3);

% Show cropped version for normalized fish 
% (only for plotting for now, because drawing on GUI relies on the true
% coordinates)

if isRefAnat && isPopout,
    % crop lengthwise (most of tail), and scale k_zres
    y_range = 81:1000;%81:990;% 81:1104; % max 1406
    z_range = z_range_ventral:345; % max 138*k_zres
    dimv_yx3 = size(anat_YX(y_range,:,:));
    dimv_yz3 = [length(y_range),length(z_range),3];
    dimv_zx3 = [length(z_range),size(anat_YX,2),3];
    
    dim_totimage = [dimv_yx3(1)+dimv_zx3(1)+10,dimv_yx3(2)+dimv_yz3(2)+10,3];
else
    dim_totimage = [dimv_yx(1)+dimv_zx(1)*k_zres_ratio+10,dimv_yx(2)+dimv_yz(2)*k_zres_ratio+10,3];    
end
tot_image=ones(dim_totimage);

% create masks (~linearized indices), to draw cells on image (later)
circlemaskIX = MakeCircularMask(radius_xy,dimv_yx);
yzmaskIX = MakeSquareMask(width_z,thickness_z,dimv_yz);
zxmaskIX = MakeSquareMask(thickness_z,width_z,dimv_zx);

% set transparancy
alpha_max = 0.4;%0.3-min(length(cIX)/1000/100,0.1);
if isWeighAlpha,    
    wIX = getappdata(hfig,'wIX');
    alpha = wIX*alpha_max;
else
    alpha = ones(size(cIX))*alpha_max;
end

%% Draw each cell on map (all 3 views)
% main loop
cIX_abs = absIX(cIX);
M_xyz = CellXYZ(cIX_abs,:);

anat_YX = DrawMasksInRGB(anat_YX,M_xyz(:,[1,2]),circlemaskIX,clrmap,gIX,alpha);
anat_YZ = DrawMasksInRGB(anat_YZ,M_xyz(:,[1,3]),yzmaskIX,clrmap,gIX,alpha/2);
anat_ZX = DrawMasksInRGB(anat_ZX,M_xyz(:,[3,2]),zxmaskIX,clrmap,gIX,alpha/2);

%% Draw masks on anat projections
if isShowMasks,  
    nMasks = length(Msk_IDs);
    % set params
    rng(1);
%     cmap2 = rand(nMasks,3);
    cmap2 = hsv(nMasks);%jet(nMasks);
    outline_radius = 3;
    msk_alpha = 0.1;%0.15 %0.25 for all projections
    white_alpha = 0.15;%0.25 % 0 for all projections
    for i = 1:nMasks,
        clr = cmap2(i,:);
        % get masks from database
        if Msk_IDs(i)==0,
            msk = full(newMask);
        else
            msk = full(MASKs.MaskDatabase(:,Msk_IDs(i)));
        end
        mask_3D = reshape(msk, [MASKs.height, MASKs.width, MASKs.Zs]);
        
        % Y-X
        masks_XY = max(mask_3D,[],3);
        masks_fit = logical(imresize(masks_XY,[size(anat_yx,1),size(anat_yx,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), strel('disk',outline_radius));
        end
        anat_YX = DrawMasksInRGB(anat_YX,[],masks_fit,clr,1,msk_alpha,white_alpha);

        % Y-Z
        masks_YZ = squeeze(max(mask_3D,[],2));
        masks_fit = logical(imresize(masks_YZ,[size(anat_yz,1),size(anat_yz,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), ones(5,2));%strel('line',outline_radius,90));
        end
        anat_YZ = DrawMasksInRGB(anat_YZ,[],masks_fit,clr,1,msk_alpha,white_alpha);

        % Z-X
        masks_ZX = squeeze(max(mask_3D,[],1))';  % notice the transpose
        masks_fit = logical(imresize(masks_ZX,[size(anat_zx,1),size(anat_zx,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), ones(2,5));%strel('line',outline_radius,0));
        end
        anat_ZX = DrawMasksInRGB(anat_ZX,[],masks_fit,clr,1,msk_alpha,white_alpha);
    end
end

%% rescale (low-res) z dimension
for k=1:3
    anat_yz2(:,:,k) = imresize(anat_YZ(:,:,k), [dimv_yz(1), dimv_yz(2)*k_zres_ratio],'nearest');   
    anat_zx2(:,:,k) = imresize(anat_ZX(:,:,k), [dimv_zx(1)*k_zres_ratio, dimv_zx(2)],'nearest');
end

% Draw scale bar
if isRefAnat
    if ~isPopout,
        y_range = 1:size(anat_YX,1);
    end
    k_um = 0.798; % ZBrain: x/y/z = 0.798/0.798/2um
    scalebar_len = round(50/k_um);
    scalebar_x_range = (600-scalebar_len+1):600;
    scalebar_y = y_range(end)-30-y_range(1)+1;
    scalebar_y_range = (scalebar_y-3):(scalebar_y+3);
    
    temp_YX = anat_YX(y_range,:,:);
    temp_YX(scalebar_y_range,scalebar_x_range,:) = 1;
else
    k_um = 0.406; % xy_res
    scalebar_len = round(50/k_um);
    
    x_end = size(anat_YX,2)-30;
    scalebar_x_range = (x_end-scalebar_len+1):x_end;
    scalebar_y = size(anat_YX,1)-30+1;
    scalebar_y_range = (scalebar_y-3):(scalebar_y+3);
    
    temp_YX = anat_YX;
    temp_YX(scalebar_y_range,scalebar_x_range,:) = 1;
end

% compile 3 projections to image panel
if isRefAnat && isPopout,    
    tot_image(dimv_zx3(1)+11:end,1:dimv_yx3(2),:) = temp_YX;%anat_YX(y_range,:,:);
    tot_image(dimv_zx3(1)+11:end,dimv_yx3(2)+11:end,:) = anat_yz2(y_range,z_range,:);
    tot_image(1:dimv_zx3(1),1:dimv_zx3(2),:) = flipud(anat_zx2(z_range,:,:));    
else
    tot_image(dimv_zx(1)*k_zres_ratio+11:end,1:dimv_yx(2),:) = temp_YX;
    tot_image(dimv_zx(1)*k_zres_ratio+11:end,dimv_yx(2)+11:end,:) = anat_yz2;
    tot_image(1:dimv_zx(1)*k_zres_ratio,1:dimv_zx(2),:) = flipud(anat_zx2);
end

% OPTIONAL BRIGHTENING
% tot_image = tot_image*2;


tot_image(tot_image(:)>1) = 1;
tot_image(tot_image(:)<0) = 0;

image(tot_image);
axis image;axis off

end