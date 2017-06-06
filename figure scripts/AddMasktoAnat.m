
folder = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\figures\Fig2\tiff stack';
tiffdir = fullfile(folder,'PT_SwimLR_allfish.tiff');
%
isPlotfromtiffstack = 1;
if isPlotfromtiffstack
    n = length(imfinfo(tiffdir));
    IM_full = cell(1,n);
    for i = 1:n
        im = double(imread(tiffdir,i))./255;
        IM_full{i} = im;%(317:1236,1:621,:);
    end
end

i_reg = 1;
range_im = 1:n;%,5:7];%[1:3,5:18];
cellarray = IM_full(i_reg,range_im);

%%
k_scale = 0.8; % this changes for every reg pair...
k_contrast = 1;

[h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
% imwrite(im_avr, fullfile(outputDir,'White_1reg_allfish_avr.tiff'), 'compression','none','writemode','overwrite');

%%
% im_avr = PT_SwimLR_avr;
im_avr = im_avr(1:1216,:,:);
%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% draw anat masks
opts.isShowMasks = 1;
opts.isShowMskOutline = 1;
opts.Msk_IDs = [232,193,194];

cIX = [];
gIX = [];
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,[],[],opts);
[h,im_mask] = DrawCellsOnAnat(I);

%% max projection
IM = cat(4,im_avr,im_mask);
mip = max(IM, [], 4);
figure;
imagesc(mip)


