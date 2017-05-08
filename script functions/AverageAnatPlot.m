function [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale)
n = length(cellarray);

sum = cellarray{1};
for i = 2:n
    im = cellarray{i};
    if isempty(im)
        disp(['image #',num2str(i),' is empty!']);
        h_anat = [];
        im_avr = [];
        return
    end
    sum = sum+im.^k_contrast;
end
im_avr = sum/sqrt(n)*k_scale;

h_anat = figure('Position',[600,50,458.5,608],'color',[1 1 1]);
set(gca,'Position',[0.01, 0.01, 0.98, 0.98]);

imagesc(im_avr)
axis equal
axis off
    
disp(['n = ' num2str(n)])
end