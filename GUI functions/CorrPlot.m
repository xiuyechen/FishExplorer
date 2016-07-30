function CorrPlot(coeffs,isPlotText,ylabels)
im = coeffs;

% red-white-blue colormap
cmap = zeros(64,3);
cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
minlim = -1; %min(min(im));
maxlim = 1; %max(max(im));

RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

image(RGB); axis equal; axis tight;
set(gca,'XTick',1:length(im),'YTick',1:length(im));

if exist('isPlotText','var'),
    if isPlotText,
        for i = 1:size(im,1),
            for j = 1:size(im,2),
                text(i-0.3, j, num2str(round(im(j,i)*100)/100));%, 'Units', 'data') % flipped i/j necessary
            end
        end
    end
end
if exist('ylabels','var'),
    set(gca,'YTickLabel',ylabels);
end
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end