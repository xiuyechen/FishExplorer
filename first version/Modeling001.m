function Modeling001(hfig)
cIX_0 = getappdata(hfig,'cIX');
gIX_0 = getappdata(hfig,'gIX');
numU = getappdata(hfig,'numU');
M = getappdata(hfig,'M');
C = FindCentroid(gIX_0,M);
[~,~,I_stimlock] = RankByStimLock_Direct(hfig,cIX_0,gIX_0,M,numU);

%% target response clusters:
% IX_rsp = [45,50,53,57,74,75]; % indices of response clusters choosen
% manually from GUI stim-lock

IX_rsp = [105,118];

%% const. params
corr_range1 = [0.5,0.9];
corr_range2 = [0.3,0.9];
thres_active = 0; % for zscored C

%% Find candidate clusters
Cz = zscore(C,0,2);

hfig_model = figure; hold on;
    
numU = length(IX_rsp); % change for combined plot...
        
for i_count = 1:length(IX_rsp),
    disp(num2str(i_count));
    i_rsp = IX_rsp(i_count);
    
    rsp = Cz(i_rsp,:);
    
    % should be lower in stim-lock ranking
    x_weights = find(I_stimlock==i_rsp);
    IX_sl = I_stimlock(1:x_weights-1);
    
    %% Method 1: filtering full centroids
    % correlation should fall within range;
    coeffs = corr(Cz(IX_sl,:)',rsp')';
    IX_cand1 = IX_sl(find(coeffs>corr_range1(1) & coeffs<corr_range1(2)));
    
    %% Method 2: masking for active periods first
    % find active periods??
    mask = rsp>thres_active;
    
    C_active = Cz(:,mask);
    rsp_active = rsp(mask);
    
    coeffs = corr(C_active(IX_sl,:)',rsp_active')';
    IX_cand2 = IX_sl(find(coeffs>corr_range2(1) & coeffs<corr_range2(2)));
    
    
    %% GLM
    IX_cand = union(IX_cand1,IX_cand2,'stable');
    
    X = Cz(IX_cand,:)';
    y = rsp';
    
    [b,dev,stats] = glmfit(X,y,'normal');
    
    yfit = glmval(b,X,'identity');
    
    numX = size(X,2);
    c = corr(y,yfit);
    
    x_weights = reshape(b(2:end),numX,[]);
    
    if 0,
        figure('Position',[50,300,1500,500]);
        
        subplot(3,1,1);
        imagesc(X');colormap gray
        
        subplot(3,1,2);
        imagesc(vertcat(y',yfit'));colormap gray
        title(['correlation = ' num2str(round(c*100)/100)]);
        
        subplot(3,1,3);hold on;
        bar(x_weights,'facecolor',[0.7 0.7 0.7]);
        xlim([1,numX])
    end
    
    %% corr-plot
    if 0,
        im = corr(horzcat(X,y,yfit));%corr(C(1,:)',C(2,:)')
        
        % red-white-blue colormap
        cmap = zeros(64,3);
        cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
        cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
        cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
        minlim = -1; %min(min(im));
        maxlim = 1; %max(max(im));
        
        RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap
        
        figure('Position',[1000,200,500,500]);
        image(RGB); axis equal; axis tight;
        
        for i = 1:size(im,2), % horizontal.. because of image axis
            for j = 1:size(im,1),
                text(i-0.3, j, num2str(round(im(i,j)*100)/100));%, 'Units', 'data')
            end
        end
    end
    %% Anatomical connection plot
    anat_yx = getappdata(hfig,'anat_yx'); % just to get size
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    CInfo = getappdata(hfig,'CInfo');
    
    range = [i_rsp; IX_cand]';
    tempI = [];
    for i = range,
        tempI = [tempI;find(gIX_0==i)];
    end
    cIX = cIX_0(tempI);
    gIX = gIX_0(tempI);
    
    % stable 'Squeeze'
    U = unique(gIX,'stable');
    numK = length(U);
    for i = 1:numK,
        old = U(i);
        gIX(gIX==old) = i;
    end

    figure(hfig_model);
    if i_count == 1,
        tot_image = DrawClustersOnMap_Modeling(CInfo,cIX,gIX,numU,i_count,anat_yx,anat_yz,anat_zx,'hsv',[1;abs(x_weights)]); 
    else
        tot_image = DrawClustersOnMap_Modeling(CInfo,cIX,gIX,numU,i_count,anat_yx,anat_yz,anat_zx,'hsv',[1;abs(x_weights)],tot_image); 
    end
    
end

end

function [C,D] = FindCentroid(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end

function [gIX,rankscore,I] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU)
periods = getappdata(hfig,'periods');

if length(periods)==1,
    period = periods{1};
    [C,~] = FindCentroid(gIX,M);
else % i_fish>=8,
    tlists = getappdata(hfig,'tlists');
    CRZt = getappdata(hfig,'CRZt');
    IX = tlists{6}; % ptomr_circ
    M_ = CRZt(cIX,IX);
    [C,~] = FindCentroid(gIX,M_);
    periods = getappdata(hfig,'periods');
    period = periods{1}+periods{2};
end
C2 = zscore(C,0,2);
C_3D = reshape(C2,size(C2,1),period,[]);
H = nanmean(nanstd(C_3D,0,3),2);
[gIX,rankscore,I] = SortH(H,gIX,numU);
end

function [gIX,B,I] = SortH(H,gIX,numU,descend) % new gIX is sorted based on H, size(H)=[numU,1];
if exist('descend','var'),
    [B,I] = sort(H,'descend');
else
    [B,I] = sort(H);
end
gIX_last = gIX;
for i = 1:numU,
    gIX(gIX_last==I(i)) = i;
end
end

function  [tot_image, dim_totimage] = DrawClustersOnMap_Modeling(CInfo,cIX,gIX,numK,i_count,anat_yx,anat_yz,anat_zx,clrmap,x_weights,tot_image_last)
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
% if ~exist('full','var'),
    displaymax = 8000;
    if length(cIX) > displaymax,
        skip = round(length(cIX)/displaymax);
        cIX = cIX(1:skip:end,:);
        gIX = gIX(1:skip:end,:);
    end
% end

if strcmp(clrmap,'jet'),
    temp = flipud(jet(numK));
else % 'hsv'
    temp = hsv(round(numK*1.1));
end
cmap = temp(1:numK,:); % extend colormap to include black

dimv_yx = size(anat_yx);
dimv_yz = size(anat_yz);
dimv_zx = size(anat_zx);
anat_YX = zeros(dimv_yx);
anat_YZ = zeros(dimv_yz);
anat_ZX = zeros(dimv_zx);
anat_ZX = flipud(anat_ZX);

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
w_white = 0.3;

for j=1:length(cIX)
    if ~isempty(CInfo(cIX(j)).center),
        %% Y-X
        %     cinds = sub2ind([dim_y,dim_x],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).center(2));
        cinds=(CInfo(cIX(j)).center(2)-1)*dimv_yx(1)+CInfo(cIX(j)).center(1); % faster equivalent, lin px idx
        labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yx(1)*dimv_yx(2)); % within bounds
        ix = gIX(j);%find(U==gIX(j));U = 1:numK;
        ixs = cinds+circle_inds(labelinds);
        anat_YX(ixs) = cmap(i_count,1)*weight*x_weights(ix) + anat_YX(ixs)*(1-weight) + w_white*(ix==1); % R
        ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2);
        anat_YX(ixs) = cmap(i_count,2)*weight*x_weights(ix) + anat_YX(ixs)*(1-weight) + w_white*(ix==1); % G
        ixs = cinds+circle_inds(labelinds)+dimv_yx(1)*dimv_yx(2)*2;
        anat_YX(ixs) = cmap(i_count,3)*weight*x_weights(ix) + anat_YX(ixs)*(1-weight) + w_white*(ix==1); % B
        
        % Y-Z
        zweight = weight/2;
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        cinds=(CInfo(cIX(j)).slice-1)*dimv_yz(1)+CInfo(cIX(j)).center(1);
        labelinds=find((cinds+yzplane_inds)>0 & (cinds+yzplane_inds)<=dimv_yz(1)*dimv_yz(2));
        ixs = cinds+yzplane_inds(labelinds);
        anat_YZ(ixs) = cmap(i_count,1)*zweight*x_weights(ix) + anat_YZ(ixs)*(1-zweight) + w_white*(ix==1);
        ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2);
        anat_YZ(ixs) = cmap(i_count,2)*zweight*x_weights(ix) + anat_YZ(ixs)*(1-zweight) + w_white*(ix==1);
        ixs = cinds+yzplane_inds(labelinds)+dimv_yz(1)*dimv_yz(2)*2;
        anat_YZ(ixs) = cmap(i_count,3)*zweight*x_weights(ix) + anat_YZ(ixs)*(1-zweight) + w_white*(ix==1);
        
        % Z-X
        zweight = weight/2;
        %     cinds = sub2ind([dim_y,dim_z],cell_info(cellsIX(j)).center(1),cell_info(cellsIX(j)).slice);
        cinds=(CInfo(cIX(j)).center(2)-1)*dimv_zx(1) +(CInfo(cIX(j)).slice);
        labelinds=find((cinds+zxplane_inds)>0 & (cinds+zxplane_inds)<=dimv_zx(1)*dimv_zx(2));
        ixs = cinds+zxplane_inds(labelinds);
        anat_ZX(ixs) = cmap(i_count,1)*zweight*x_weights(ix) + anat_ZX(ixs)*(1-zweight) + w_white*(ix==1);
        ixs = cinds+zxplane_inds(labelinds)+dimv_zx(1)*dimv_zx(2);
        anat_ZX(ixs) = cmap(i_count,2)*zweight*x_weights(ix) + anat_ZX(ixs)*(1-zweight) + w_white*(ix==1);
        ixs = cinds+zxplane_inds(labelinds)+dimv_zx(1)*dimv_zx(2)*2;
        anat_ZX(ixs) = cmap(i_count,3)*zweight*x_weights(ix) + anat_ZX(ixs)*(1-zweight) + w_white*(ix==1);
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

if exist('tot_image_last','var'),
    tot_image = tot_image + tot_image_last;
    tot_image(tot_image>1) = 1;
end

image(tot_image);
axis image;axis ij;axis off

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

% function [gIX, numK] = SqueezeGroupIX(gIX)
% U = unique(gIX);
% numK = length(U);
% for i = 1:numK,
%     old = U(i);
%     gIX(gIX==old) = i;
% end
% end