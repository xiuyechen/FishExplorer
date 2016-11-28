function I = LoadMultiFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap)
currdir = GetCurrentDataDir();
load(fullfile(currdir,'BasicMultiFishInfo.mat'));


I = [];

I.isRefAnat = 1;%getappdata(hfig,'isRefAnat');
I.isPopout = getappdata(hfig,'isPopout');

cIX_m = GetMultifishCellIndex(cIX,i_fish);

offset = length(cell2mat(CellXYZ_multi(1:i_fish-1)));

I.CellXYZ = getappdata(hfig,'CellXYZ_norm');
I.anat_yx = getappdata(hfig,'anat_yx_norm');
I.anat_yz = getappdata(hfig,'anat_yz_norm');
I.anat_zx = getappdata(hfig,'anat_zx_norm');

I.absIX = getappdata(hfig,'absIX');
if exist('cIX_plot','var'),
    I.cIX = cIX_plot;
else
    I.cIX = getappdata(hfig,'cIX');
    cIX = I.cIX;
end
if exist('gIX_plot','var'),
    I.gIX = gIX_plot;
else
    I.gIX = getappdata(hfig,'gIX');
    gIX = I.gIX;
end

% get colormap
if exist('clrmap','var'),
    I.clrmap = clrmap;
else
    % get colormap
    clrmap_name = getappdata(hfig,'clrmap_name');
    numK = getappdata(hfig,'numK');
    numK = double(max(numK,max(gIX)));% sadly this is not always true
    I.clrmap = GetColormap(clrmap_name,numK);
end

% anat masks
I.isShowMasks = getappdata(hfig,'isShowMasks');
if I.isShowMasks,
    isShowMskOutline = getappdata(hfig,'isShowMskOutline');
    MASKs = getappdata(hfig,'MASKs');
    Msk_IDs = getappdata(hfig,'Msk_IDs');
    I.Msk_IDs = Msk_IDs;
    if Msk_IDs == 0,
        I.newMask = getappdata(hfig,'newMask');
    end
end

% set transparancy
isWeighAlpha = getappdata(hfig,'isWeighAlpha');
alpha_max = 0.4;%0.3-min(length(cIX)/1000/100,0.1);
if isWeighAlpha,
    wIX = getappdata(hfig,'wIX');
    I.clr_alpha = wIX*alpha_max;
else
    I.clr_alpha = ones(size(cIX))*alpha_max;
end

% fish outline
I.isShowFishOutline = getappdata(hfig,'isShowFishOutline');
I.FishOutline = getappdata(hfig,'FishOutline');
    if I.isShowFishOutline,
        FishOutline = getappdata(hfig,'FishOutline');
        % option with some background showing
%         anat_yx = 0.5*anat_yx + 0.5*repmat(FishOutline.outline_XY,[1,1,3]);
%         anat_yz = 0.5*anat_yz + 0.5*repmat(FishOutline.outline_YZ,[1,1,3]);
%         anat_zx = 0.5*anat_zx + 0.5*repmat(FishOutline.outline_ZX,[1,1,3]);
        
        I.anat_yx = repmat(FishOutline.outline_XY,[1,1,3]);
        I.anat_yz = repmat(FishOutline.outline_YZ,[1,1,3]);
        I.anat_zx = repmat(FishOutline.outline_ZX,[1,1,3]);
    end
    
end


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