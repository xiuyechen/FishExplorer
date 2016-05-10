function  [tot_image, dim_totimage] = DrawCellsOnAnatProj(hfig,isRefAnat,isPopout)
%% load
if ~isRefAnat,
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres = 20;
    radius_xy = 7;
    width_z = 10;
    thickness_z = 1;
else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_yx = getappdata(hfig,'anat_yx_norm');
    anat_yz = getappdata(hfig,'anat_yz_norm');
    anat_zx = getappdata(hfig,'anat_zx_norm');
    k_zres = 2.5;
    radius_xy = 3;
    width_z = 5;
    thickness_z = 3;
end

absIX = getappdata(hfig,'absIX');
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
numK = getappdata(hfig,'numK');

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

% down-sample
if ~isPopout,
    displaymax = 8000;
    if length(cIX) > displaymax,
        skip = round(length(cIX)/displaymax);
        cIX = cIX(1:skip:end,:);
        gIX = gIX(1:skip:end,:);
    end
end

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

% create masks (~linearized indices), to draw cells on image (later)
circlemaskIX = MakeCircularMask(radius_xy,dimv_yx);
yzmaskIX = MakeSquareMask(width_z,thickness_z,dimv_yz);
zxmaskIX = MakeSquareMask(thickness_z,width_z,dimv_zx);

% set transparancy
alpha_max = 0.3-min(length(cIX)/1000/100,0.1);
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

anat_YX = DrawMasksInRGB(anat_YX,M_xyz(:,[1,2]),circlemaskIX,cmap,gIX,alpha);
anat_YZ = DrawMasksInRGB(anat_YZ,M_xyz(:,[1,3]),yzmaskIX,cmap,gIX,alpha/2);
anat_ZX = DrawMasksInRGB(anat_ZX,M_xyz(:,[3,2]),zxmaskIX,cmap,gIX,alpha/2);

%% Draw masks on anat projections
if isShowMasks,  
    nMasks = length(Msk_IDs);
    % set params
    cmap2 = jet(nMasks);
    outline_radius = 3;
    msk_alpha = 0.25;
    white_alpha = 0.05;
    M_xyz = [];
    for i = 1:nMasks,
        clr = cmap2(i,:);
        % get masks from database
        msk = full(MASKs.MaskDatabase(:,Msk_IDs(i)));
        mask_3D = reshape(msk, [MASKs.height, MASKs.width, MASKs.Zs]);
        
        % Y-X
        masks_XY = max(mask_3D,[],3);
        masks_fit = logical(imresize(masks_XY,[size(anat_yx,1),size(anat_yx,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), strel('disk',outline_radius));
        end
        anat_YX = DrawMasksInRGB(anat_YX,M_xyz,masks_fit,clr,1,msk_alpha,white_alpha);

        % Y-Z
        outline_radius = 6;
        masks_YZ = squeeze(max(mask_3D,[],2));
        masks_fit = logical(imresize(masks_YZ,[size(anat_yz,1),size(anat_yz,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), ones(5,2));%strel('line',outline_radius,90));
        end
        anat_YZ = DrawMasksInRGB(anat_YZ,M_xyz,masks_fit,clr,1,msk_alpha,white_alpha);

        % Z-X
        masks_ZX = squeeze(max(mask_3D,[],1))';  % notice the transpose
        masks_fit = logical(imresize(masks_ZX,[size(anat_zx,1),size(anat_zx,2)]));
        if isShowMskOutline,
            masks_fit = imdilate(edge(masks_fit), ones(2,5));%strel('line',outline_radius,0));
        end
        anat_ZX = DrawMasksInRGB(anat_ZX,M_xyz,masks_fit,clr,1,msk_alpha,white_alpha);
    end
end

%% rescale (low-res) z dimension
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