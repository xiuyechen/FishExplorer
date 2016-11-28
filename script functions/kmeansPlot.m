function kmeansPlot(M,groupIX,Centroids)
% rng('default');
% [groupIX,C] = kmeans(M,numK);
% [groupIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);
if exist('Centroids','var'),
    groupIX = HierClus_Direct(Centroids,groupIX);
end
numK = length(unique(groupIX));

% sort data matrix based on clustering
[~,I] = sort(groupIX);
im = M;
im = im(I,:);

%% convert imagesc effect into rgb matrix 'RGB'
cmap = gray(64);
im = mat2gray(im);
minlim = 0;
maxlim = 1;
RGB = ImageToRGB(im,cmap,minlim,maxlim);

%% add vertical color code
nLines = size(M,1);
barwidth = max(round(size(M,2)/30),1);

% make colormap
clrmap = hsv(round(numK*1.1)); % 'stretch out' hsv -> less circular

% find group divisions for color-bar
idx = groupIX(I);
ix_div = [find(diff(idx));length(idx)];

bars = ones(nLines,barwidth);
if numK>1,
    bars(1:ix_div(1),:) = idx(ix_div(1));
    for i = 2:length(ix_div),
        % paint color
        bars(ix_div(i-1)+1:ix_div(i),:) = idx(ix_div(i));
    end
end
im_bars = reshape(clrmap(bars,:),[size(bars),3]);

% add a white division margin
div = ones(nLines,round(barwidth/2),3);

% put 3 parts together
im = horzcat(RGB,div,im_bars);

%% plot cluster division lines
figure;
imagesc(im); hold on
axis off
% label axes
% if isObjectSpace,
%     xlabels = [Bstar.bodyPart,{''},{''}]; % adjust for the colorbar and margin
%     set(gca,'XTickLabel',xlabels,'XTickLabelRotation',45);
%     set(gca,'XTick',1:(size(im,2)+2));
%     set(gca,'TickLength',[0 0]);
%     
%     ylabel('Objects');
%     set(gca,'YTick',[]); 
% else
%     set(gca,'YTick',1:size(im,1));
%     set(gca,'YTickLabel',Bstar.bodyPart(I));
%     set(gca,'TickLength',[0 0]);
%     
%     xlabel('Objects');
% end

% plot cluster division lines in red
if numK>1,
    for i = 1:length(ix_div),% = numK-1,
        y = ix_div(i)+0.5;
        plot([0.5,size(im,2)+0.5],[y,y],'r','Linewidth',0.5);
    end
end
