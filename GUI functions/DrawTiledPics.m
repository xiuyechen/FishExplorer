function DrawTiledPics(hfig)
absIX = getappdata(hfig,'absIX');
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
anat_stack = getappdata(hfig,'anat_stack');

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

isShowMasks = getappdata(hfig,'isShowMasks');
if isShowMasks,
    MASKs = getappdata(hfig,'MASKs');
    Msk_IDs = getappdata(hfig,'Msk_IDs');    
end

% timestamp = datestr(now,'mmddyy_HHMMSS');
% tiffName = ['stack_' timestamp '.tif'];

U = unique(gIX);
numK = length(U);

anat_stack2=zeros(size(anat_stack,1), size(anat_stack,2), size(anat_stack,3) ,3);
nPlanes=size(anat_stack,3);
dimv_yxz=size(anat_stack);
stacklen=numel(anat_stack);

% circle=makeDisk2(10,21);
radius = 10; dim = 21;
center=floor(dim/2)+1;
circle=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<=radius
            circle(x,y)=1;
        end
    end
end


[r, v]=find(circle);
r=r-11;v=v-11;
circle_inds  = r*dimv_yxz(1)+v;
cmap = hsv(round(numK*1.1));
% cmap = [0.3 1 0];
weight = 0.3;

for i=1:nPlanes,
    anat_stack2(:,:,i,:)=repmat(imNormalize99(anat_stack(:,:,i))/4,[1 1 1 3]);
end

% for j=1:length(cIX)
%     cinds=(CInfo(cIX(j)).center(2)-1)*dimv_yxz(1)+CInfo(cIX(j)).center(1);
%     labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yxz(1)*dimv_yxz(2));
%     zinds=dimv_yxz(1)*dimv_yxz(2)*(CInfo(cIX(j)).slice-1);
%     ix = find(U==gIX(j));
%     ixs = cinds+circle_inds(labelinds)+zinds;
%     anat_stack2(ixs)=cmap(ix,1)*weight + anat_stack2(ixs)*(1-weight);
%     ixs = cinds+circle_inds(labelinds)+zinds+stacklen;
%     anat_stack2(ixs)=cmap(ix,2)*weight + anat_stack2(ixs)*(1-weight);
%     ixs = cinds+circle_inds(labelinds)+zinds+stacklen*2;
%     anat_stack2(ixs)=cmap(ix,3)*weight + anat_stack2(ixs)*(1-weight);
% end

    cIX_abs = absIX(cIX);
for j=1:length(cIX)        
    if ~isempty(CellXYZ(cIX_abs(j),1)),
        ix = gIX(j);
    %% Y-X        
        cinds=(CellXYZ(cIX_abs(j),2)-1)*dimv_yx(1)+CellXYZ(cIX_abs(j),1); % linear pixel index, faster equivalent of:
        %     cinds = sub2ind([dim_y,dim_x],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).center(2));
        labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yx(1)*dimv_yx(2)); % within bounds
        ixs = cinds+circle_inds(labelinds);
        anat_YX(ixs) = cmap(ix,1)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % R
        ixs = cinds+circle_inds(labelinds) + dimv_yx(1)*dimv_yx(2);
        anat_YX(ixs) = cmap(ix,2)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % G
        ixs = cinds+circle_inds(labelinds) + dimv_yx(1)*dimv_yx(2)*2;
        anat_YX(ixs) = cmap(ix,3)*alpha(j) + anat_YX(ixs)*(1-alpha(j)); % B
    end
end

if isShowMasks,  
    nMasks = length(Msk_IDs);
    % set params
    cmap2 = jet(nMasks);
    msk_alpha = 0.3;
    white_alpha = 0.1;
    
    for i = 1:nMasks,
        msk = full(MASKs.MaskDatabase(:,Msk_IDs(i)));
        mask_3D = reshape(msk, [MASKs.height, MASKs.width, MASKs.Zs]);
        
        % Y-X
%         masks_XY = max(mask_3D,[],3);
%         masks_fit = logical(imresize(masks_XY,[size(anat_yx,1),size(anat_yx,2)]));
        ixs = find(masks_fit);
        anat_stack2(ixs) = (cmap2(i,1)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % R;
        ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2);
        anat_stack2(ixs) = (cmap2(i,2)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % G;
        ixs = find(masks_fit) + dimv_yx(1)*dimv_yx(2)*2;
        anat_stack2(ixs) = (cmap2(i,3)*msk_alpha + anat_stack2(masks_fit)*(1-msk_alpha))*(1-white_alpha) + white_alpha; % B;
    end
end
%% view stack sequentially
% figure;
% for i_plane = 1:nPlanes,
%     im = squeeze(ave_stack2(:,:,i_plane,:));
%     image(im);
%     axis image; axis off
%     pause(0.1)
% end

% lay out

h = figure('Position',[10,20,round(1136*1.4),round(683*1.4)],'color','k');
hold on;

axes('Position',[0,0,1,1]); % BS fix to get the figure background to stay when saving pics
image(zeros(3,3,3));axis off

numRow = 3;
numP = nPlanes - mod(nPlanes,numRow);
numCol = numP/numRow;

for i_plane = 1:numP,
    im = squeeze(anat_stack2(:,:,i_plane,:));
    im = imresize(im,0.25);
    
    [col, row] = ind2sub([numCol,numRow],i_plane); % this is intensionally inverted...
    posVec = [(col-1)*1/numCol,1-row*1/numRow,1/numCol,1/numRow];
    axes('Position',posVec);
    image(im);
    axis image; axis off
    %     pause(0.1)
end

end