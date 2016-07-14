M_dir = GetFishDirectories();
save_masterdir = GetNestedDataDir();
 
range_fish = GetFishRange();
%%
range_fish = 6:15;
figure
 
for i_fish = range_fish,
     data_dir = M_dir{i_fish};
        
    %% init
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(save_dir,'CoreInfo.mat'));
    load(fullfile(save_dir,'OptionalInfo.mat'));
    
    %% get anatomy projections
    % x-y view
    im = max(anat_stack,[],3);
    out=imNormalize99(im);
    anat_yx = repmat(out,[1 1 3]);
    
    % y-z view
    im = squeeze(max(anat_stack,[],2));
    out=imNormalize99(im);
    anat_yz = repmat(out,[1 1 3]);
    
    % x-z view
    im = squeeze(max(anat_stack,[],1));
    out=imNormalize99(im);
    out = flipud(out'); %%%% empirically necessary...
    anat_zx = repmat(out,[1 1 3]);

   
    
    %% test Plot: all antomy outliers
    subplot(1,2,1)
    title(num2str(i_fish));
    cIX = IX_inval;
    gIX = (1:length(cIX))';
    numK = round(cIX(end)/1000)+1;
    BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);
    
    %% test Plot: all remaining cells
    subplot(1,2,2)
    I_v_Holder = ones(1,numcell_full);
    I_v_Holder(IX_inval) = 0;
    cIX = find(I_v_Holder);
    gIX = (1:length(cIX))';
    numK = round(cIX(end)/1000)+1;
    BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);
    
    waitforbuttonpress
end