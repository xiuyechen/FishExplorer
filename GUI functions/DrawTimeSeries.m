function DrawTimeSeries(hfig,cIX_plot,gIX_plot,h1,isPopout,isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS)
if ~exist('cIX_plot','var'),
    cIX = getappdata(hfig,'cIX');
else
    cIX = cIX_plot;
end
if ~exist('gIX_plot','var'),
    gIX = getappdata(hfig,'gIX');
else
    gIX = gIX_plot;
end
if ~exist('h1','var'),
    h1 = gca;
end
if isempty(h1),
    h1 = gca;
end
if ~exist('isPopout','var'),
    isPopout = getappdata(hfig,'isPopout');
end
if ~exist('isCentroid','var'),
    isCentroid = getappdata(hfig,'isCentroid');
end
if ~exist('isPlotLines','var'),
    isPlotLines = getappdata(hfig,'isPlotLines');
end
if ~exist('isPlotBehavior','var'),
    isPlotBehavior = getappdata(hfig,'isPlotBehavior');
end
if ~exist('isPlotRegWithTS','var'),
    isPlotRegWithTS = getappdata(hfig,'isPlotRegWithTS');
end

% load
numK = getappdata(hfig,'numK');
if isempty(numK),
    numK = max(gIX);
end
behavior = getappdata(hfig,'behavior');
stim = getappdata(hfig,'stim');
clrmap_name = getappdata(hfig,'clrmap_name');
rankscore = getappdata(hfig,'rankscore');
rankID = getappdata(hfig,'rankID');
fpsec = getappdata(hfig,'fpsec');
iswrite = (rankID>=2);

nCells = length(cIX);
regressor = GetRegressor(hfig);
i_fish = getappdata(hfig,'i_fish');

isPlotRegSameRow = 0;
if ~exist('isPlotRegWithTS','var'),
    isPlotRegWithTS = 0;
end    

M = getappdata(hfig,'M');

nFrames = size(M,2);

if isPopout, % ...for isCentroid, skipping before getting Centroid causes trouble...
%     down-sample
    ds_cap = 1000;
    nCells = length(cIX);
    skip = max(1,round(nCells/ds_cap));
    M = M(1:skip:end,:); % down-sized
    gIX = gIX(1:skip:end,:); % down-sized
end

axis off;
pos = get(h1,'Position'); % [left, bottom, width, height]    

% double-check numK
numK = max(max(gIX),numK);

if isPlotLines,
    if 1, % (just to collapse isPlotLines = true)               
        %% settings
        len = pos(3);
%         fpsec = getappdata(hfig,'fpsec');%1.97;
        
        % set position grid
        nClus = length(unique(gIX));
        
        if isPlotBehavior,
            nExtraRows = 3; % number of extra rows
        else
            nExtraRows = 2;
        end        
        if isPlotRegWithTS && ~isPlotRegSameRow,
            nLines = nClus + 1;
            nRows = max(8,nLines)+nExtraRows+1;
        else
            nLines = nClus;
            nRows = max(8,nLines)+nExtraRows;
        end
        lineH = pos(4)/nRows;
        xv = 1:nFrames;        
                
        %% Set colormap        
        clrmap = GetColormap(clrmap_name,numK)*0.8; % make darker by *0.8
        
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
        U = unique(gIX);
        
        j_pos = 0;
        if isPlotRegWithTS,
            if ~isPlotRegSameRow,
                j_pos = 1;
            end
            
            j = 1;
            % draw regressor
            sub_pos = [pos(1),pos(2)+pos(4)-lineH*(j+1),len,lineH*0.95];
            axes('Position',sub_pos); hold on;
            plot(xv,regressor,'-','Linewidth',1,'color','k')
            axis tight;axis off
        end
        
        for j = 1:nClus,
            k = U(j);
            Ys = M(gIX==k,:);            
            ymean = mean(Ys,1);
            ySTD=std(Ys,0,1); % /sqrt(size(Ys,1))??
            ySTD_upper=ymean+ySTD;
            ySTD_lower=ymean-ySTD;
            
            % Draw
            sub_pos = [pos(1),pos(2)+pos(4)-lineH*(j+1+j_pos),len,lineH*0.95];
            axes('Position',sub_pos); hold on;
            % draw std
            h = fill([xv fliplr(xv)], [ySTD_upper fliplr(ySTD_lower)],0.7+0.2*clrmap(k,:));
            set(h, 'EdgeColor','none')
            % draw mean
            plot(xv,ymean,'-','Linewidth',1,'color',clrmap(k,:))
            axis tight;axis off
        end
        
        %% plot scale bar
        axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+2),len,lineH]); hold on;
        if nFrames<600, % <~10 min, plot 1min scale bar
            plot([1,10*fpsec],[0.75,0.75],'k-');
            text(pos(1),0.5,'10 sec','HorizontalAlignment','left')
        elseif nFrames<2400, % <~10 min, plot 1min scale bar
            plot([1,60*fpsec],[0.75,0.75],'k-');
            text(pos(1)+60*fpsec/2,0.5,'1 min','HorizontalAlignment','center')
        else % >~10 min, plot 5 min scale bar
            plot([1,300*fpsec],[0.75,0.75],'k-');
            text(pos(1)+300*fpsec/2,0.5,'5 min','HorizontalAlignment','center')
        end
        xlim([xv(1),xv(end)]);ylim([0,1]);
        axis off
                
        %% draw behavior bar
        if isPlotBehavior,
            h = axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+3),len,0.7/nRows]);    
            DrawBehaviorBar(h,behavior,i_fish);
            
            axes('Position',[0.02,pos(2)+pos(4)-lineH*(nLines+3),pos(1)-0.02,0.7/nRows]);
            DrawArrowsIcon(isPopout);
        end

    end
    
else % ~isPlotLines, i.e. plot all traces as grayscale map
    if isCentroid,
        M_ = FindClustermeans(gIX,M);
        gIX_ = unique(gIX);
    else
        M_ = M;
        gIX_ = gIX;
    end
    
    %% set size of stimulus-bar relative to whole plot
    % each stim bar (if applicable) takes up 1 unit, and behavior bar takes up 2.
    if ~isPopout,
        barratio = 0.025;
    else
        barratio = 0.05;
    end
    
    %% Prepare cluster data   
    % sort traces by index
    im = M_;
    [nLines,nFrames] = size(im);
    [~,I] = sort(gIX_);
    im = im(I,:);
    
    %% convert imagesc effect into rgb matrix 'RGB'
    cmap = gray(64);
    im = AutoScaleImage0to1(im); % scaling to min/max 0/1
    minlim = 0;
    maxlim = 1;
    
    RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap
    
    %% add vertical color code
     
    clrmap0 = GetColormap(clrmap_name,numK);
    clrmap = vertcat(clrmap0,[0,0,0]); % extend colormap to include black
    bwidth = max(round(nFrames/30),1);
    
    idx = gIX_(I);
    ix_div = [find(diff(idx));length(idx)];
    
    bars = ones(nLines,bwidth);
    if numK>1,
        bars(1:ix_div(1),:) = idx(ix_div(1));
        for i = 2:length(ix_div),
            % paint color
            bars(ix_div(i-1)+1:ix_div(i),:) = idx(ix_div(i));
        end
    end
    im_bars = reshape(clrmap(bars,:),[size(bars),3]);
    % add a white division margin
    div = ones(nLines,round(bwidth/2),3);
    % put 3 parts together
    im = horzcat(RGB,div,im_bars);
    
    %% plot main figure    
    [s1,s2,s3] = size(im);
       
    % plot 'im'
    hold off;
    if isPopout,
        if ~isPlotBehavior,
            set(h1,'Position',[pos(1),pos(2),pos(3),pos(4)*(1-barratio)]);
        else
            set(h1,'Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*(1-barratio*3)]);
        end
    else % in main window, isPlotBehavior is always true
        % but need to leave extra space for the bottom axis (compared to isPopout)
        set(h1,'Position',[pos(1),pos(2)+pos(4)*barratio*3,pos(3),pos(4)*(1-barratio*4)]);
    end
    
    % (if too few rows, pad to 30 rows with white background)
    if s1<30,
        temp = im;
        im = ones(30,s2,s3);
        im(1:s1,:,:) = temp;
    end    
    image(im);
    colormap gray;
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
    
    % label axes
    if ~exist('nCells','var'),
        nCells = nLines;
    end
    ylabel(['Number of cells: ' num2str(nCells)]);
    set(gca,'YTick',[],'XTick',[]);
    set(gcf,'color',[1 1 1]);
            
    % write in text labels (numbers corresponding to vertical colorbar)
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
    roughhalfbarheight = 100;
    [stimbar_,halfbarheight] = GetStimBar(roughhalfbarheight,stim);
    
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
    if ~isPopout, % plot stimbar in the bottom too
        axes('Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*barratio]);
        image(stimbar);axis off;hold on;
        plot([1,length(stim)],[size(stimbar,1),size(stimbar,1)],'k','Linewidth',0.5); % plot bottom border
        plot([0.5,0.5],[0,size(stimbar,1)],'k','Linewidth',0.5); % plot left border
        plot([length(stim),length(stim)],[0,size(stimbar,1)],'k','Linewidth',0.5); % plot right border
    end
    
    %% Draw behavior
    if isPlotBehavior,
        h = axes('Position',[pos(1),pos(2),pos(3),pos(4)*barratio*2]);        
        DrawBehaviorBar(h,behavior,i_fish,barlength);        
        axes('Position',[0.015,pos(2),0.03,pos(4)*barratio*2]);
        DrawArrowsIcon(isPopout);
    else
        h = h1;
    end

%     % label horizontal axis: approximate time length
%     set(h,'XColor',[1 1 1]);
%     if s2/2/60<5, % <5 min
%         xlabel(['Time ~ ' num2str(round(s2/fpsec)) 'sec'],'Color',[0 0 0]);
%     else
%         xlabel(['Time ~ ' num2str(round(s2/fpsec/60)) 'min'],'Color',[0 0 0]);
%     end
    
    %% Draw fish icon:
    if ~isPopout,
        axes('Position',[0.015, pos(2)+pos(4)*barratio*2, 0.03, pos(4)*barratio]);
    else
        axes('Position',[0.015,pos(2)+pos(4)*(1-barratio),0.03,pos(4)*barratio]);
    end
    DrawFishIcon;
end

set(gcf,'Color','w');
end

%% Internal Functions

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

function DrawBehaviorBar(h,behavior,i_fish,barlength)
regressors = GetMotorRegressor(behavior,i_fish);
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
if 1, % plot all 5 lines
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
    m = vertcat(temp(1,:),temp(2,:),temp(3,:));    
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
