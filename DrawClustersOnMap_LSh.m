function  [tot_image, dim_totimage] = DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,clrmap,full)
%% formatting
[s1,s2] = size(cIX);
if s2>s1,
    cIX = cIX';
end
[s1,s2] = size(gIX);
if s2>s1,
    gIX = gIX';
end

% down-sample
if ~exist('full','var'),
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

if strcmp(clrmap,'jet'),
    temp = flipud(jet(numK));
else % 'hsv'
    temp = hsv(round(numK*1.1));
end
cmap = temp(1:numK,:); % extend colormap to include black

anat_YX = anat_yx/4;
anat_YZ = anat_yz/4;
anat_ZX = anat_zx/4;
anat_ZX = flipud(anat_ZX);
dimv_yx = size(anat_YX);
dimv_yz = size(anat_YZ);
dimv_zx = size(anat_ZX);

k_zres = 20;
anat_yz2=zeros(dimv_yz(1),dimv_yz(2)*k_zres,3);
anat_zx2=zeros(dimv_zx(1)*k_zres,dimv_zx(2),3);
dim_totimage = [dimv_yx(1)+dimv_zx(1)*k_zres+10,dimv_yx(2)+dimv_yz(2)*k_zres+10,3];
tot_image=ones(dim_totimage);
% tot_image=zeros(dim_totimage);
% tot_image(:,dimv_yx(2)+(1:10),:)=1;

% find index manipulation vector to darw circle
circle=makeDisk2(7,15); % make mask of filled circle % (7,15)
mask = zeros(dimv_yx(1),dimv_yx(2));
mask(1:15,1:15) = circle;
ix = find(mask);
cix = sub2ind([dimv_yx(1),dimv_yx(2)],8,8);% 8
circle_inds = ix - cix;

yzplane_inds = -5:5;
zxplane_inds = -5*dimv_zx(1):dimv_zx(1):5*dimv_zx(1);

weight = 0.3-min(length(cIX)/1000/100,0.1);
for j=1:length(cIX)
    if ~isempty(CInfo(cIX(j)).center),
        %% Y-X
        %     cinds = sub2ind([dim_y,dim_x],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).center(2));
        cinds=(CInfo(cIX(j)).center(2)-1)*dimv_yx(1)+CInfo(cIX(j)).center(1); % faster equivalent, lin px idx
        labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yx(1)*dimv_yx(2)); % within bounds
        ix = gIX(j);%find(U==gIX(j));U = 1:numK;
        ixs = cinds+circle_inds(labelinds);
        anat_YX(ixs) = cmap(ix,1)*weight + anat_YX(ixs)*(1-weight); % R
        ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2);
        anat_YX(ixs) = cmap(ix,2)*weight + anat_YX(ixs)*(1-weight); % G
        ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2)*2;
        anat_YX(ixs) = cmap(ix,3)*weight + anat_YX(ixs)*(1-weight); % B
        
        % Y-Z
        zweight = weight/2;
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        cinds=(CInfo(cIX(j)).slice-1)*dimv_yz(1)+CInfo(cIX(j)).center(1);
        labelinds=find((cinds+yzplane_inds)>0 & (cinds+yzplane_inds)<=dimv_yz(1)*dimv_yz(2));
        ixs = cinds+yzplane_inds(labelinds);
        anat_YZ(ixs) = cmap(ix,1)*zweight + anat_YZ(ixs)*(1-zweight);
        ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2);
        anat_YZ(ixs) = cmap(ix,2)*zweight + anat_YZ(ixs)*(1-zweight);
        ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2)*2;
        anat_YZ(ixs) = cmap(ix,3)*zweight + anat_YZ(ixs)*(1-zweight);
                
        % Z-X
        zweight = weight/2;
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        cinds=(CInfo(cIX(j)).center(2)-1)*dimv_zx(1) +(CInfo(cIX(j)).slice);
        labelinds=find((cinds+zxplane_inds)>0 & (cinds+zxplane_inds)<=dimv_zx(1)*dimv_zx(2));
        ixs = cinds+zxplane_inds(labelinds);
        anat_ZX(ixs) = cmap(ix,1)*zweight + anat_ZX(ixs)*(1-zweight);
        ixs = cinds+zxplane_inds(labelinds)+dimv_zx(1)*dimv_zx(2);
        anat_ZX(ixs) = cmap(ix,2)*zweight + anat_ZX(ixs)*(1-zweight);
        ixs = cinds+zxplane_inds(labelinds)+dimv_zx(1)*dimv_zx(2)*2;
        anat_ZX(ixs) = cmap(ix,3)*zweight + anat_ZX(ixs)*(1-zweight);
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