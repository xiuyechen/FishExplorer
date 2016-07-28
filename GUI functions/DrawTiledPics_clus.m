function DrawTiledPics_clus(hfig)
% isRefAnat = getappdata(hfig,'isRefAnat');
% if ~isRefAnat,
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_stack = getappdata(hfig,'anat_stack');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres = 20;
    radius_xy = 7;
    width_z = 10;
    thickness_z = 1;
% else
%     CellXYZ = getappdata(hfig,'CellXYZ_norm');
%     anat_stack = getappdata(hfig,'anat_stack_norm');
%     anat_yx = getappdata(hfig,'anat_yx_norm');
%     anat_yz = getappdata(hfig,'anat_yz_norm');
%     anat_zx = getappdata(hfig,'anat_zx_norm');
%     k_zres = 2.5;
%     radius_xy = 3;
%     width_z = 5;
%     thickness_z = 3;
% end

% isShowMasks = getappdata(hfig,'isShowMasks');
% if isShowMasks,
%     MASKs = getappdata(hfig,'MASKs');
%     Msk_IDs = getappdata(hfig,'Msk_IDs');    
% end

% load params
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
absIX = getappdata(hfig,'absIX');

segment = 8;
% nPlanes = size(anat_stack,3);
nPlanes = ceil(max(gIX)/segment);

% initialize image stack
dimv_yxz = size(anat_stack);
anat_stack2 = zeros([dimv_yxz(1:2),nPlanes,3]);

% make cell-shaped circular mask
radius_xy = 7;
circlemaskIX = MakeCircularMask(radius_xy,dimv_yxz(1:2));

% make color-map
numK = segment; %length(U);
clrmap = hsv(round(double(numK)*1.1));
% set transparency
alpha = ones(size(cIX))*0.5;

%% main: coloring of stack
cIX_abs = absIX(cIX);
M_xyz = CellXYZ(cIX_abs,:);
stack_alpha = 0.25;
for i = 1:nPlanes,
    anat_plane = stack_alpha*repmat(squeeze(imNormalize99(anat_stack(:,:,i))),[1 1 3]);
    IX = find(ceil(gIX/segment)==i); % for DrawTiledPics.m: IX = find(M_xyz(:,3)==i);
    clrIX = gIX(IX)-(i-1)*segment;
    if ~isempty(IX),
        anat_stack2(:,:,i,:) = DrawMasksInRGB(anat_plane,M_xyz(IX,[1,2]),circlemaskIX,clrmap,clrIX,alpha(IX));
    end
end

%% lay out

h = figure('Position',[10,20,round(1136*1.4),round(683*1.4)],'color','k'); 
hold on;

axes('Position',[0,0,1,1]); % BS fix to get the figure background to stay when saving pics
image(zeros(3,3,3));axis off

numRow = 3;
numP = ceil(max(gIX)/segment); %nPlanes - mod(nPlanes,numRow);
numCol = ceil(numP/numRow);

for i_plane = 1:numP,
    im = squeeze(anat_stack2(:,:,i_plane,:));
%     im = imresize(im,0.25);
    
    [col, row] = ind2sub([numCol,numRow],i_plane); % this is intensionally inverted...
    posVec = [(col-1)*1/numCol,1-row*1/numRow,1/numCol,1/numRow];
    axes('Position',posVec);
    image(im);
    axis image; axis off
%     pause(0.1)
end
end
