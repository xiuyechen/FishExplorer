function DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite)
pos = get(gca,'Position');
barratio = 0.03;

%% Prepare cluster data
% down-sample
displaymax = 1000;
numcell = size(M,1);
if numcell > displaymax,
    skip = round(numcell/displaymax);
    M = M(1:skip:end,:);
    gIX = gIX(1:skip:end,:);
end

numK = double(max(double(numK),double(max(gIX))));

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
plot([0.5,s2+0.5],[0.5,0.5],'k','Linewidth',0.5);
if numK>1,
    for i = 1:length(ix_div),% = numK-1,
        y = ix_div(i)+0.5;
        plot([0.5,s2+0.5],[y,y],'k','Linewidth',0.5);
    end
end

ylabel(['Number of cells: ' num2str(numcell)]);
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

%% Draw stim bars
halfbarheight = 100;
stimbar_ = GetStimBar(halfbarheight,stim);

% pad with margin on right to fit combo-image size
barlength = size(im,2);
stimbar = ones(halfbarheight*2,barlength,3);
stimbar(:,1:nFrames,:) = stimbar_(:,1:nFrames,:); % crop to size

% top axis
axes('Position',[pos(1),pos(2)+pos(4)*(1-barratio),pos(3),pos(4)*barratio]);
image(stimbar);axis off;hold on;
plot([1,length(stim)],[1,1],'k','Linewidth',0.5); % plot top border
plot([0.5,0.5],[0,size(stimbar,1)],'k','Linewidth',0.5); % plot left/right border
plot([length(stim),length(stim)],[0,size(stimbar,1)],'k','Linewidth',0.5);

% bottom axis
axes('Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*barratio]);
image(stimbar);axis off;hold on;
plot([1,length(stim)],[size(stimbar,1),size(stimbar,1)],'k','Linewidth',0.5); % plot bottom border
plot([0.5,0.5],[0,size(stimbar,1)],'k','Linewidth',0.5); % plot left/right border
plot([length(stim),length(stim)],[0,size(stimbar,1)],'k','Linewidth',0.5);

%% Draw fictive
axes('Position',[pos(1),pos(2),pos(3),pos(4)*barratio*2]);

temp = ones(size(fictive,1),barlength);
temp(:,1:length(fictive)) = fictive;
imagesc(temp);colormap hot
set(gca,'YTick',[],'XTick',[]);
set(gcf,'color',[1 1 1]);
set(gca, 'box', 'off')
hold on;axis ij;

% plot division lines
for i = 0:3,
    y = i+0.5;
    plot([0.5,length(fictive)+0.5],[y,y],'w','Linewidth',0.5);
end
% labels
names = {'Left','Right','Forward','Raw L','Raw R'};
x = -s2*0.05;
for i = 1:5,
    y = i;
    text(x,y,names{i},'Fontsize',7);
end

%% horizontal axis: approximate time length, assuming 2 fps (real is 1.97 for now)
set(gca,'XColor',[1 1 1]);
if s2/2/60<5, % <5 min
    l = xlabel(['Length of display: ~ ' num2str(s2/2) 'sec'],'Color',[0 0 0]);
else
    l = xlabel(['Length of display: ~ ' num2str(s2/2/60) 'min'],'Color',[0 0 0]);
end

%% Draw fish icon:
if exist('fishpic.jpg','file')==2,
    fishpic = imread('fishpic.jpg');
    axes('Position',[0.015, 0.095, 0.03, 0.02]);
    image(fishpic); axis off; axis equal
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
