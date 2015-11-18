function BasicPlotMaps(cIX,gIX,M,CellXYZ,photostate,anat_yx,anat_yz,anat_zx) %(cIX,gIX,cIX_0,M_0,CIF,numK,M)
% M = M_0(cIX_0(cIX),:);
numK = length(unique(gIX));
%%
dataFR = 1;
figure('Position',[100 50 1350 900]);%,'DeleteFcn',@closefigure_Callback);
h1 = axes('Position',[0.05, 0.06, 0.5, 0.85]); % left ~subplot
BasicDrawClusters(h1,M,gIX,dataFR,numK,photostate);

h2 = axes('Position',[0.58, 0.06, 0.40, 0.85]); % right ~subplot
BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);%,anat_zx,'hsv');
% set(gcf,'PaperUnits', 'inches', 'PaperSize', [12, 6])

end

function BasicDrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite)
pos = get(gca,'Position');
barratio = 0.02;
numK = double(max(double(numK),double(max(gIX))));

% down-sample
displaymax = 1000;
numcell = size(M,1);
if numcell > displaymax,
    skip = round(numcell/displaymax);
    M = M(1:skip:end,:);
    gIX = gIX(1:skip:end,:);
end

% sort traces by index
im = M;
[nLines,nFrames] = size(im);
[~,I] = sort(gIX);
im = im(I,:);

%% convert imagesc effect into rgb matrix 'RGB'
if dataFR, % is drawing cell-response fluorescence
    cmap = gray(64);
    im = AutoScaleImage0to1(im); % scaling to min/max 0/1
    minlim = 0;
    maxlim = 1;
else % create blue-white-red map, for regression results
    cmap = zeros(64,3);
    cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
    cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
    cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
    minlim = -1; %min(min(im));
    maxlim = 1; %max(max(im));
end
RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

%% add vertical color code
if exist('clrmap','var'),
    if strcmp(clrmap,'jet'),
        temp = flipud(jet(numK));
    else % 'hsv'
        temp = hsv(round(numK*1.1));
    end
else % 'hsv'
    temp = hsv(round(numK*1.1));
end
cmap2 = vertcat(temp(1:numK,:),[0,0,0]); % extend colormap to include black
bwidth = max(round(nFrames/30),1);

idx = gIX(I);
ix_div = [find(diff(idx));length(idx)];

bars = ones(nLines,bwidth);
if numK>1,
    bars(1:ix_div(1),:) = idx(ix_div(1));
    for i = 2:length(ix_div),
        % paint color
        bars(ix_div(i-1)+1:ix_div(i),:) = idx(ix_div(i));
    end
end
im_bars = reshape(cmap2(bars,:),[size(bars),3]);
% add a white division margin
div = ones(nLines,round(bwidth/2),3);
% put 3 parts together
im = horzcat(RGB,div,im_bars);

%% plot figure
% h_fig = figure('Position',[100 0 1300 900],'Name',['round ' num2str(i_step)]); %[left, bottom, width, height]
[s1,s2,s3] = size(im);

hold off;
set(h1,'Position',[pos(1),pos(2)+pos(4)*barratio*3,pos(3),pos(4)*(1-4*barratio)]); 
if s1<30,
    temp = im;
    im = ones(30,s2,s3);
    im(1:s1,:,:) = temp;
end
image(im);
set(gca, 'box', 'off')

hold on;

% plot cluster division lines
if numK>1,
    for i = 1:length(ix_div),% = numK-1,
        y = ix_div(i)+0.5;
        plot([0.5,s2+0.5],[y,y],'k','Linewidth',0.5);
%         plot([1,s2*kx],[y,y],'k','Linewidth',0.5);
    end
end

ylabel(['ROI number (downsampled): ' num2str(s1)]);
set(gca,'YTick',[],'XTick',[]);
set(gcf,'color',[1 1 1]);

colormap gray;

% write in number label
x = s2*1.003;
y_last = 0;
for i = 1:length(ix_div),
    % avoid label crowding
    margin = 0.015*s1;
    y0 = ix_div(i)+0.5;
    y = max(y_last+margin,y0);
    if i<length(ix_div),
        ynext = ix_div(i+1);
    else
        ynext = y0+margin*2;
    end
    % draw if not squeezed too much away
    if y < y0+margin*2 && y < 0.5*(y0+ynext), %nLines+stimheight-margin,
        if exist('iswrite','var') && iswrite,
            try
                text(x,y,[num2str(idx(ix_div(i))) ': ' num2str(rankscore(i))],'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
                    'FontUnits','normalized','FontSize',0.015,'Tag','label');
            catch
                disp('rankscore invalid');
                text(x,y,num2str(idx(ix_div(i))),'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
                    'FontUnits','normalized','FontSize',0.015,'Tag','label');
            end
            
        else
            text(x,y,num2str(idx(ix_div(i))),'HorizontalAlignment','Left','VerticalAlignment','Bottom',...
                'FontUnits','normalized','FontSize',0.015,'Tag','label');
        end
        y_last = y;
    end
end

if exist('stim','var'),
    barlength = size(im,2);
    DrawBars(pos,barratio,nFrames,barlength,stim);
    
    if exist('fictive','var'),
        % Draw fictive
        axes('Position',[pos(1),pos(2),pos(3),pos(4)*barratio*2]);
        temp = ones(size(fictive,1),barlength);
        temp(:,1:length(fictive)) = fictive;
        imagesc(temp);colormap hot
        axis off;
    end
end

end

%% Internal Functions

function im = AutoScaleImage0to1(im)
temp = mat2gray(im);
low_high = stretchlim(temp,[0.005,0.995]);
temp(temp<low_high(1)) = low_high(1);
temp(temp>low_high(2)) = low_high(2);
im = mat2gray(temp);
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end

function DrawBars(pos,barratio,nFrames,barlength,photostate) % horizontal stimulus bar
% pos = get(gca,'Position'); % [0.05, 0.06, 0.5, 0.85]

% global htop hbottom;
barheight = 100;
stimheight = 2*barheight;

m1 = 0.3*ones(barheight,length(photostate));
m2 = m1; % bottom half
x = photostate;
    
%     2,3,4,5,12,13,14,15,99
%     3,1,3,2,4 ,10,11,12,
    
if max(x)<=3,
    % 0 = all black; 1 = black/white; 2 = white/black; 3 = all white; (4 = all gray;)
    % 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
    
    % top (~ projection right)
    m1(:,x==0 | x==2 | x==5) = 0; % black
    m1(:,x==4 | x==6 | x==7) = 0.5; % grey
    m1(:,x==1 | x==3 | x==8) = 1; % white
    % bottom (~ projection left)
    m2(:,x==0 | x==1 | x==7) = 0;
    m2(:,x==4 | x==5 | x==8) = 0.5;
    m2(:,x==2 | x==3 | x==6) = 1;
    
    temp = repmat(vertcat(m1,m2),[1 1 3]); % 3 color layers in dim3
else
    %%
    % 2 = white/white
    % 3 = black/white
    % 4 = white/white
    % 5 = white/black
    %
    % 12 = gray
    % 13 = forward grating (very slow, more for calibration)
    % 14 = rightward grating
    % 15 = leftward grating  (* I hope I got right/left correct - if you think it's flipped it probably is)
    %
    % 99 = spontaneous
    % top (~ projection right)
    
    m1(:,x==3) = 0; % black
    m1(:,x==2 | x==4 | x==5) = 1; % white
    % bottom (~ projection left)
    m2(:,x==5) = 0;
    m2(:,x==2 | x==3 | x==4) = 1;
    
    temp = repmat(vertcat(m1,m2),[1 1 3]); % 3 color layers in dim3
    %%
    temp(:,x==14,1) = 1; % red
    temp(:,x==13,2) = 1; % green
    temp(:,x==15,3) = 1; % blue
    
    % grey
    temp(:,x==12,1) = 0.8;
    temp(:,x==12,2) = 0.8;
    temp(:,x==12,3) = 0.8;
    
    % yellow
    temp(:,x==99,1) = 1;
    temp(:,x==99,2) = 1;
    temp(:,x==99,3) = 0.6;
end

n = ceil(nFrames/length(photostate)); % tile, to size or bigger
temp = repmat(temp,1,n);
stimbar = ones(stimheight,barlength,3);
stimbar(:,1:nFrames,:) = temp(:,1:nFrames,:); % crop to size

% if isempty(htop),
axes('Position',[pos(1),pos(2)+pos(4)*(1-barratio),pos(3),pos(4)*barratio]);
% elseif ~ishandle(htop),
%     htop = axes('Position',[pos(1),pos(2)+pos(4)*(1-barratio),pos(3),pos(4)*barratio]);
% else axes(htop)
% end
image(stimbar);axis off;hold on;
% plot top border
plot([1,length(photostate)],[1,1],'k','Linewidth',0.5);

% if isempty(hbottom),
axes('Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*barratio]);
% elseif ~ishandle(hbottom),
%     hbottom = axes('Position',[pos(1),pos(2),pos(3),pos(4)*barratio]);
% else axes(hbottom)
% end

image(stimbar);axis off;hold on;
% plot bottom border
plot([1,length(photostate)],[stimheight,stimheight],'k','Linewidth',0.5);
%%
end
