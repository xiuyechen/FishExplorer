function WriteZstack(hfig,tiffdir,cIX,gIX,clrmap)
if ~exist('cIX','var')
    cIX = getappdata(hfig,'cIX');
end
if ~exist('gIX','var')
    gIX = getappdata(hfig,'gIX');
end
if ~exist('clrmap','var')
    numK = double(max(gIX));
    clrmap_name = getappdata(hfig,'clrmap_name');
    clrmap = GetColormap(clrmap_name,numK); %hsv(round(double(numK)*1.1));
end

if ~exist('tiffdir','var')
    % get save path
    timestamp = datestr(now,'mmddyy_HHMMSS');
    tiffName = ['stack_' timestamp '.tif'];
    [file,path] = uiputfile(tiffName,'Save tif stack');
    tiffdir = fullfile(path,file);
end

% load params
absIX = getappdata(hfig,'absIX');
isRefAnat = getappdata(hfig,'isRefAnat');
if ~isRefAnat
    CellXYZ = getappdata(hfig,'CellXYZ');    
    anat_stack = getappdata(hfig,'anat_stack');
    radius_xy = 7;
else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_stack = getappdata(hfig,'anat_stack_norm');
    radius_xy = 4;
end

% initialize image stack
anat_stack2 = zeros([size(anat_stack),3]);
nPlanes = size(anat_stack,3);
dimv_yxz = size(anat_stack);

% make cell-shaped circular mask
circlemaskIX = MakeCircularMask(radius_xy,dimv_yxz(1:2));

% make color-map

% set transparency
alpha = ones(size(cIX))*0.5;

% main: coloring of stack
cIX_abs = absIX(cIX);
M_xyz = CellXYZ(cIX_abs,:);
stack_alpha = 0.25;
for i = 1:nPlanes,
    anat_plane = stack_alpha*repmat(imNormalize99(anat_stack(:,:,i)),[1 1 1 3]);
    IX = find(M_xyz(:,3)==i);
    if ~isempty(IX),
        anat_stack2(:,:,i,:) = DrawMasksInRGB(anat_plane,M_xyz(IX,[1,2]),circlemaskIX,clrmap,gIX(IX),alpha(IX));
    end
end

% display each plane and save as tif
h = figure;
for i_plane = 1:nPlanes,
    im = squeeze(anat_stack2(:,:,i_plane,:));
    image(im);axis equal; axis off
    drawnow;
    % save tiff
    if (i_plane == 1)
        imwrite(im, tiffdir, 'compression','none','writemode','overwrite')
    else
        imwrite(im, tiffdir, 'compression','none','writemode','append')
    end
    %     pause(0.2)
end
close(h)

end