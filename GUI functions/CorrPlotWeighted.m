function cmap = CorrPlotWeighted(coeffs,gIX)%,isPlotText,ylabels)

coeffs(isnan(coeffs)) = 0;

% size of 'squares' corresponds to cluster size (from 'sizes')
% im = coeffs;
[gIX, numK] = SqueezeGroupIX(gIX);
U = unique(gIX);
sizes = zeros(1,numK);
for i = 1:numK
    IX = find(gIX==U(i));
    sizes(i) = length(IX);
end

nCells = length(gIX);

if nCells<500
    im = zeros(nCells,nCells);
    for i = 1:numK
        for j = 1:numK
            IX1 = find(gIX==U(i));
            IX2 = find(gIX==U(j));
            [X,Y] = ndgrid(IX1,IX2);
            im(X,Y) = coeffs(i,j);
        end
    end
else % downsample
    skip = floor(nCells/500);
    gIX_ds = gIX(1:skip:end);
    nCells_ds = length(gIX_ds);
    im = zeros(nCells_ds,nCells_ds);
    
    for i = 1:numK
        for j = 1:numK
            IX1 = find(gIX_ds==U(i));
            IX2 = find(gIX_ds==U(j));
            [X,Y] = ndgrid(IX1,IX2);
            im(X,Y) = coeffs(i,j);
        end
    end
end
       
% red-white-blue colormap
cmap = zeros(64,3);
cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
minlim = -1; %min(min(im));
maxlim = 1; %max(max(im));

RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

image(RGB); axis equal; axis tight;
if length(im)<50
    set(gca,'XTick',1:length(im),'YTick',1:length(im));
else
    set(gca,'XTick',[],'YTick',[]);
end
set(gca,'TickLength',[0,0]);

% if exist('isPlotText','var'),
%     if isPlotText,
%         for i = 1:size(im,1),
%             for j = 1:size(im,2),
%                 text(i-0.3, j, num2str(round(im(j,i)*100)/100));%, 'Units', 'data') % flipped i/j necessary
%             end
%         end
%     end
% end
% if exist('ylabels','var'),
%     set(gca,'YTickLabel',ylabels);
% end
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end