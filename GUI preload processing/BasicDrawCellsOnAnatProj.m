function  [tot_image, dim_totimage] = BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz)
%% formatting
[s1,s2] = size(cIX);
if s2>s1,
    cIX = cIX';
end
[s1,s2] = size(gIX);
if s2>s1,
    gIX = gIX';
end

% get numK
if exist('numK','var'),
    numK = double(max(numK,max(gIX)));
else
    numK = double(max(gIX));
end

temp = hsv(round(numK*1.1));
cmap = temp(1:numK,:); % extend colormap to include black

anat_YX = anat_yx;
anat_YZ = anat_yz;
dimv_yx = size(anat_YX);
dimv_yz = size(anat_YZ);

k_zres = 20;
anat_yz2=zeros(dimv_yz(1),dimv_yz(2)*k_zres,3);
% anat_yz2_ori=anat_yz2;
dim_totimage = [dimv_yx(1),dimv_yx(2)+dimv_yz(2)*k_zres+10,3];
tot_image=zeros(dim_totimage);
tot_image(:,dimv_yx(2)+(1:10),:)=1;

% find index manipulation vector to darw circle
circle=makeDisk2(7,15); % make mask of filled circle % (7,15)
mask = zeros(dimv_yx(1),dimv_yx(2));
mask(1:15,1:15) = circle;
ix = find(mask);
cix = sub2ind([dimv_yx(1),dimv_yx(2)],8,8);% 8
circle_inds = ix - cix;

yzplane_inds = -5:5;

weight = 0.5-min(length(cIX)/1000/100,0.3);
for j=1:length(cIX)
    %     cinds = sub2ind([dim_y,dim_x],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).center(2));
    cinds=(CellXYZ(cIX(j),2)-1)*dimv_yx(1)+CellXYZ(cIX(j),1); % faster equivalent, lin px idx
    labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yx(1)*dimv_yx(2)); % within bounds
    ix = gIX(j);%find(U==gIX(j));U = 1:numK;
    ixs = cinds+circle_inds(labelinds);
    anat_YX(ixs) = cmap(ix,1)*weight + anat_YX(ixs)*(1-weight); % R
    ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2);
    anat_YX(ixs) = cmap(ix,2)*weight + anat_YX(ixs)*(1-weight); % G
    ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2)*2;
    anat_YX(ixs) = cmap(ix,3)*weight + anat_YX(ixs)*(1-weight); % B
    
    zweight = weight/2;
%     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
    cinds=(CellXYZ(cIX(j),3)-1)*dimv_yz(1)+CellXYZ(cIX(j),1);
    labelinds=find((cinds+yzplane_inds)>0 & (cinds+yzplane_inds)<=dimv_yz(1)*dimv_yz(2));
    ixs = cinds+yzplane_inds(labelinds);
    anat_YZ(ixs) = cmap(ix,1)*zweight + anat_YZ(ixs)*(1-zweight);
    ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2);
    anat_YZ(ixs) = cmap(ix,2)*zweight + anat_YZ(ixs)*(1-zweight);
    ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2)*2;
    anat_YZ(ixs) = cmap(ix,3)*zweight + anat_YZ(ixs)*(1-zweight);
end
% rescale (low-res) z dimension
for k=1:3
    anat_yz2(:,:,1)=imresize(anat_YZ(:,:,1), [dimv_yz(1) dimv_yz(2)*k_zres],'nearest');
    anat_yz2(:,:,2)=imresize(anat_YZ(:,:,2), [dimv_yz(1) dimv_yz(2)*k_zres],'nearest');
    anat_yz2(:,:,3)=imresize(anat_YZ(:,:,3), [dimv_yz(1) dimv_yz(2)*k_zres],'nearest');
end

tot_image(:,1:dimv_yx(2),:)=anat_YX;
tot_image(:,dimv_yx(2)+11:end,:)=anat_yz2;

tot_image(tot_image(:)>1)=1;
tot_image(tot_image(:)<0)=0;

image(tot_image);
axis equal;axis off

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