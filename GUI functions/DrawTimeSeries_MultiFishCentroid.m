function DrawTimeSeries_MultiFishCentroid(hfig,AllCentroids,range_fish,range_clus,isPopout)
% load
numK = getappdata(hfig,'numK');
behavior = getappdata(hfig,'behavior');
stim = getappdata(hfig,'stim');
clrmap = getappdata(hfig,'clrmap');
% rankscore = getappdata(hfig,'rankscore');
rankID = getappdata(hfig,'rankID');
% iswrite = (rankID>=2);
cIX = getappdata(hfig,'cIX');
% nCells = length(cIX);

% M = getappdata(hfig,'M');
for i_fish = 1:range_fish,
C = AllCentroids.Centroids{i_fish};
gIX = getappdata(hfig,'gIX');
nFrames = size(M,2);

% if isPopout, % ...for isCentroid, skipping before getting Centroid causes trouble...
%     %     down-sample
%     ds_cap = 1000;
%     nCells = length(cIX);
%     skip = max(1,round(nCells/ds_cap));
%     M = M(1:skip:end,:); % down-sized
%     gIX = gIX(1:skip:end,:); % down-sized
% end

axis off;
pos = get(h1,'Position'); % [left, bottom, width, height]

% double-check numK
numK = max(max(gIX),numK);

% if isPlotLines,
%% settings
len = pos(3);
fpsec = getappdata(hfig,'fpsec');%1.97;

% set position grid
nLines = length(range_clus{i_fish});%length(unique(gIX));
if isPlotBehavior,
    nExtraRows = 3; % number of extra rows
else
    nExtraRows = 2;
end
nRows = max(8,nLines)+nExtraRows;
lineH = pos(4)/nRows;
xv = 1:nFrames;

%% Set colormap
if strcmp(clrmap,'jet'),
    clr = flipud(jet(double(numK)))*0.8; % make darker by *0.8
else % if strcmp(clrmap,'hsv'),
    clr = hsv(round(double(numK)*1.1))*0.7; % make darker by *0.7
end
% clr = lines(nLines);
% clr = ones(nLines,3)*0.3; % uniform charcoal

%% draw stimulus bar
stimbar = GetStimBar(1,stim);
sub_pos = [pos(1),pos(2)+pos(4)-3/4*lineH,len,lineH/2];
h = axes('Position',sub_pos);

image(stimbar);
set(h,'Xtick',[],'Ytick',[])

% draw fish icon
sub_pos = [pos(1)*0.2,pos(2)+pos(4)-3/4*lineH,pos(1)*0.7,lineH/2];
axes('Position',sub_pos);
DrawFishIcon;

%% draw cluster means
% U = unique(gIX);

for j = range_clus{i_fish},%1:nLines,
%     k = U(j);
%     Ys = M(gIX==k,:);
%     ymean = mean(Ys);
%     ySTD=std(Ys); % /sqrt(size(Ys,1))??
%     ySTD_upper=ymean+ySTD;
%     ySTD_lower=ymean-ySTD;
    ymean = C(j,:);

    % Draw
    sub_pos = [pos(1),pos(2)+pos(4)-lineH*(j+1),len,lineH];
    axes('Position',sub_pos); hold on;
    % draw std
%     h = fill([xv fliplr(xv)], [ySTD_upper fliplr(ySTD_lower)],0.5+0.3*clr(k,:));
%     set(h, 'EdgeColor','none')
    % draw mean
    plot(xv,ymean,'-','Linewidth',1,'color',clr(k,:))
    axis tight;axis off
end

%% plot scale bar
axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+2),len,lineH]); hold on;
if nFrames<1200, % <~10 min, plot 1min scale bar
    plot([1,60*fpsec],[0.75,0.75],'k-');
    text(pos(1)+60*fpsec/2,0.5,'1 min','HorizontalAlignment','center')
else % >~10 min, plot 10min scale bar
    plot([1,600*fpsec],[0.75,0.75],'k-');
    text(pos(1)+600*fpsec/2,0.5,'10 min','HorizontalAlignment','center')
end
xlim([xv(1),xv(end)]);ylim([0,1]);
axis off

%% draw behavior bar
if isPlotBehavior,
    h = axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+3),len,0.7/(nLines+nExtraRows)]);
    DrawBehaviorBar(h,behavior);
    
    axes('Position',[0.02,pos(2)+pos(4)-lineH*(nLines+3),pos(1)-0.02,0.7/(nLines+nExtraRows)]);
    DrawArrowsIcon(isPopout);
end

end

function im = AutoScaleImage0to1(im)
im(isnan(im)) = 0;
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

function DrawFishIcon
if exist('fishpic.jpg','file')==2,
    pic = imread('fishpic.jpg');
    fishpic = imresize(pic,0.5);
    
    image(fishpic); axis off; axis equal; axis tight
end
end

function DrawArrowsIcon(isPopout)
if ~isPopout,
    if exist('arrows.jpg','file')==2,
        pic = imread('arrows.jpg');
        arrowpic = imresize(pic,0.2);
        image(arrowpic); axis off; axis image; axis tight
    end
else
    if exist('arrows_wide.jpg','file')==2,
        pic = imread('arrows_wide.jpg');
        arrowpic = imresize(pic,0.2);
        image(arrowpic); axis off; axis image; axis tight
    end
end
end

function DrawBehaviorBar(h,behavior,barlength)
regressors = GetMotorRegressor(behavior);
m = reshape([regressors(1:5).im],length(behavior),[])';
m = AutoScaleImage0to1(m);

% additional line-by-line normalization! (optional)
for i = 1:5,
    temp = m(i,:);
    m(i,:) = (temp-min(temp))/(max(temp)-min(temp));
end

if exist('barlength','var'),
    temp = ones(size(m,1),barlength);
    temp(:,1:length(m)) = m;
else
    temp = m;
end
if 0, % plot all 5 lines
    imagesc(temp);colormap gray
    set(h,'YTick',[],'XTick',[]);
    %     set(gcf,'color',[1 1 1]);
    set(h, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
    for i = 0:3,
        y = i+0.5;
        plot([0.5,length(behavior)+0.5],[y,y],'w','Linewidth',0.5);
    end
    % labels
%     names = {'Left','Right','Forward','Raw L','Raw R'};
%     x = -s2*0.05;
%     for i = 1:5,
%         y = i;
%         text(x,y,names{i},'Fontsize',7);
%     end
    
else % only plot top 3 lines
    m = vertcat(temp(1,:),temp(3,:),temp(2,:));    
    imagesc(m);colormap gray
    set(h,'YTick',[],'XTick',[]);
    %     set(gcf,'color',[1 1 1]);
    set(h, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
    for i = 0:2,
        y = i+0.5;
        plot([0.5,length(behavior)+0.5],[y,y],'w','Linewidth',0.5);
    end    
end
end