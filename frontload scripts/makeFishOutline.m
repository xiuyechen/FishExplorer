% make fish outline

Msk_IDs = [275,1,94,114];
MASKs = getappdata(hfig,'MASKs');

msk = MASKs.MaskDatabase(:,Msk_IDs(1));
for i_mask = 2:length(Msk_IDs),
    msk = msk+MASKs.MaskDatabase(:,Msk_IDs(i_mask));
end

%% eyes
Msk = full(MASKs.MaskDatabase(:,78));
mask_3D_eye = reshape(Msk, [MASKs.height, MASKs.width, MASKs.Zs]);

%%
Msk = full(msk);
mask_3D = reshape(Msk, [MASKs.height, MASKs.width, MASKs.Zs]);
     
%%
mask_XY = max(mask_3D,[],3);
outline_XY_brain = imdilate(edge(logical(mask_XY)), strel('disk',3));

mask_YZ = squeeze(max(mask_3D,[],2));
outline_YZ_brain = imdilate(edge(logical(mask_YZ)), ones(5,2));%strel('line',outline_radius,90));

mask_ZX = squeeze(max(mask_3D,[],1))';  % notice the transpose
outline_ZX_brain = imdilate(edge(logical(mask_ZX)), ones(2,5));%strel('line',outline_radius,0));

%%
mask_XY_eye = max(mask_3D_eye,[],3);
outline_XY_eye = imdilate(edge(logical(mask_XY_eye)), strel('disk',1));

mask_YZ_eye = squeeze(max(mask_3D_eye,[],2));
outline_YZ_eye = imdilate(edge(logical(mask_YZ_eye)), ones(1,1));%strel('line',outline_radius,90));

mask_ZX_eye = squeeze(max(mask_3D_eye,[],1))';  % notice the transpose
outline_ZX_eye = imdilate(edge(logical(mask_ZX_eye)), ones(1,1));%strel('line',outline_radius,0));

outline_XY_eye(logical(mask_XY)) = 0; 
outline_YZ_eye(logical(mask_YZ)) = 0; 
outline_ZX_eye(logical(mask_ZX)) = 0; 

%%
outline_XY = logical(outline_XY_brain + outline_XY_eye); 
outline_YZ = logical(outline_YZ_brain + outline_YZ_eye); 
outline_ZX = logical(outline_ZX_brain + outline_ZX_eye); 

%%
% newMask = sparse(reshape(stack3_mask, [s1*s2*s3, 1]));

figure;
subplot(1,3,1)
imagesc(outline_XY);axis equal
subplot(1,3,2)
imagesc(outline_YZ);
subplot(1,3,3)
imagesc(outline_ZX);

%%
data_dir = GetCurrentDataDir();
save(fullfile(data_dir,'FishOutline.mat'),'outline_XY','outline_YZ','outline_ZX');