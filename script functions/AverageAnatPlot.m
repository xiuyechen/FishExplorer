function [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale)
n = length(cellarray);


if exist('k_contrast','var') & exist('k_scale','var')
    sum = im2double(cellarray{1});
    for i = 2:n
        im = im2double(cellarray{i});
        if isempty(im)
            disp(['image #',num2str(i),' is empty!']);
            h_anat = [];
            im_avr = [];
            return
        end
        sum = sum+im.^k_contrast;
    end
    im_avr = sum/sqrt(n)*k_scale;
else
    sum = im2double(cellarray{1});
    for i = 2:n
        im = im2double(cellarray{i});
        if isempty(im)
            disp(['image #',num2str(i),' is empty!']);
            h_anat = [];
            im_avr = [];
            return
        end
        sum = sum+im;
    end
    sum = sum/n;
    %     im_avr = imadjust(sum,[0.1,0.9],[0.1,0.9]);
    HSV = rgb2hsv(sum);
    V = squeeze(HSV(:,:,3));
    N = histcounts(V,0:0.01:1);
    thres = find(N(1:99)>100,1,'last')/100; % cutoff unit: pixels
    V2 = V/thres;
    V2(V2>1)= 1;
    HSV(:,:,3) = V2;
    im_avr = hsv2rgb(HSV);
%     im_avr = imadjust(sum,[0,thres],[0,1]);
    
end


h_anat = figure('Position',[600,50,458.5,608],'color',[1 1 1]);
set(gca,'Position',[0.01, 0.01, 0.98, 0.98]);

imagesc(im_avr)
axis equal
axis off

% disp(['n = ' num2str(n)])
text(20,1190,['n=',num2str(n)],'color','w','fontsize',15');
end
