function DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite,isPopout,isPlotLines,isPlotFictive)
axis off;
pos = get(h1,'Position'); % [left, bottom, width, height]    

if isPlotLines,
    if 1, % (just to collapse isPlotLines = true)
        %% settings
        len = pos(3);
        fps = 1.97;
        
        % set position grid
        nLines = length(unique(gIX));
        if isPlotFictive,
            nExtraRows = 3; % number of extra rows
        else
            nExtraRows = 2;
        end        
        nRows = max(8,nLines)+nExtraRows;
        lineH = pos(4)/nRows;
%         Lpos = [pos(1),pos(2)+pos(4)-lineH*nRows,pos(3),lineH*nRows];
% %         Lpos = [pos(1),pos(2)+pos(4)-lineH*(nLines+nExtraRows),pos(3),lineH*(nLines+nExtraRows)];
%         set(h1,'Position',Lpos);

        %     figure('Position',[500,900-(min(nLines,15)+nDiv)*50,600,(min(nLines,15)+nDiv)*50],'color',[1 1 1]);
        xv = 1:size(M,2);        
                
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
        rpos = [pos(1),pos(2)+pos(4)-3/4*lineH,len,lineH/2];
        h = axes('Position',rpos);
        
        image(stimbar);
        set(h,'Xtick',[],'Ytick',[])
        
        % draw fish icon
        rpos = [pos(1)*0.2,pos(2)+pos(4)-3/4*lineH,pos(1)*0.7,lineH/2];
        axes('Position',rpos);
        DrawFishIcon;    
                
        %% draw cluster means
        Ymeans = zeros(nLines,size(M,2));
        U = unique(gIX);
        for j = 1:nLines,
            k = U(j);
            %     ymean = Cz(k,:);
            Ys = M(gIX==k,:);
            ymean = mean(Ys,1);
            Ymeans(k,:) = ymean;
            rpos = [pos(1),pos(2)+pos(4)-lineH*(j+1),len,lineH];
            axes('Position',rpos); hold on;
            
            ySTD=std(Ys); % /sqrt(size(Ys,1))
            ySTD_upper=ymean+ySTD;
            ySTD_lower=ymean-ySTD;
            h = fill([xv fliplr(xv)], [ySTD_upper fliplr(ySTD_lower)],0.5+0.3*clr(k,:));
            set(h, 'EdgeColor','none')
            
            plot(xv,ymean,'-','Linewidth',1,'color',clr(k,:))
            axis tight;axis off
        end
        
        %% plot scale bar
        axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+2),len,lineH]); hold on;
        if size(M,2)<1200, % <~10 min, plot 1min scale bar
            plot([1,60*fps],[0.75,0.75],'k-');
            text(pos(1)+60*fps/2,0.5,'1 min','HorizontalAlignment','center')
        else % >~10 min, plot 10min scale bar
            plot([1,600*fps],[0.75,0.75],'k-');
            text(pos(1)+600*fps/2,0.5,'10 min','HorizontalAlignment','center')
        end
        xlim([xv(1),xv(end)]);ylim([0,1]);
        axis off
                
        %% draw fictive behavior bar
        if isPlotFictive,
            h = axes('Position',[pos(1),pos(2)+pos(4)-lineH*(nLines+3),len,0.7/(nLines+nExtraRows)]);
            DrawFictiveBar(h,fictive);
            
            axes('Position',[0.02,pos(2)+pos(4)-lineH*(nLines+3),pos(1)-0.02,0.7/(nLines+nExtraRows)]);
            DrawArrowsIcon(isPopout);
        end

    end
    
else % ~isPlotLines, i.e. plot all traces as grayscale map

    if ~isPopout,
        barratio = 0.025;
    else
        barratio = 0.05;
    end
    %% Prepare cluster data
    % down-sample
    displaymax = 1000;
    numcell = size(M,1);
    if numcell > displaymax,
        skip = round(numcell/displaymax);
        M = M(1:skip:end,:);
        gIX = gIX(1:skip:end,:);
    end
    
%     numK = double(max(double(numK),double(max(gIX))));
    
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
    if ~isPlotFictive,
        set(h1,'Position',[pos(1),pos(2),pos(3),pos(4)*(1-barratio)]);
    else
        set(h1,'Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*(1-barratio*3)]);
    end
    % axis off
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
    if ~isPopout,
        axes('Position',[pos(1),pos(2)+pos(4)*barratio*2,pos(3),pos(4)*barratio]);
        image(stimbar);axis off;hold on;
        plot([1,length(stim)],[size(stimbar,1),size(stimbar,1)],'k','Linewidth',0.5); % plot bottom border
        plot([0.5,0.5],[0,size(stimbar,1)],'k','Linewidth',0.5); % plot left/right border
        plot([length(stim),length(stim)],[0,size(stimbar,1)],'k','Linewidth',0.5);
    end
    %% Draw fictive
    if isPlotFictive,
        h = axes('Position',[pos(1),pos(2),pos(3),pos(4)*barratio*2]);        
        DrawFictiveBar(h,fictive,barlength);        
        axes('Position',[0.015,pos(2),0.03,pos(4)*barratio*2]);
        DrawArrowsIcon(isPopout);
    else
        h = h1;
    end

    % label horizontal axis: approximate time length, assuming 1.97 fps
    set(h,'XColor',[1 1 1]);
    if s2/2/60<5, % <5 min
        xlabel(['Time ~ ' num2str(round(s2/1.97)) 'sec'],'Color',[0 0 0]);
    else
        xlabel(['Time ~ ' num2str(round(s2/1.97/60)) 'min'],'Color',[0 0 0]);
    end
    
    %% Draw fish icon:
    if ~isPopout,
        axes('Position',[0.015, pos(2)+pos(4)*barratio*2, 0.03, pos(4)*barratio]);
    else
        axes('Position',[0.015,pos(2)+pos(4)*(1-barratio),0.03,pos(4)*barratio]);
    end
    DrawFishIcon;
end
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

function DrawFictiveBar(h,fictive,barlength)
regressors = GetMotorRegressor(fictive);
fc = reshape([regressors(1:5).im],length(fictive),[])';
fc = AutoScaleImage0to1(fc);
if exist('barlength','var'),
    temp = ones(size(fc,1),barlength);
    temp(:,1:length(fc)) = fc;
else
    temp = fc;
end
if 0, % plot all 5 lines
%     imagesc(temp);colormap gray
%     set(h,'YTick',[],'XTick',[]);
%     %     set(gcf,'color',[1 1 1]);
%     set(h, 'box', 'off')
%     hold on;axis ij;
%     
%     % plot division lines
%     for i = 0:3,
%         y = i+0.5;
%         plot([0.5,length(fictive)+0.5],[y,y],'w','Linewidth',0.5);
%     end
%     % labels
%     names = {'Left','Right','Forward','Raw L','Raw R'};
%     x = -s2*0.05;
%     for i = 1:5,
%         y = i;
%         text(x,y,names{i},'Fontsize',7);
%     end
    
else % only plot top 3 lines
    fc = vertcat(temp(1,:),temp(3,:),temp(2,:));
    imagesc(fc);colormap gray
    set(h,'YTick',[],'XTick',[]);
    %     set(gcf,'color',[1 1 1]);
    set(h, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
    for i = 0:2,
        y = i+0.5;
        plot([0.5,length(fictive)+0.5],[y,y],'w','Linewidth',0.5);
    end    
end
end
