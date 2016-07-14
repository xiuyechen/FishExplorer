clear all;close all;clc

M_dir = GetFishDirectories();
save_masterdir = GetNestedDataDir();

range_fish = GetFishRange();


for i_fish = range_fish,    
    disp(['i_fish = ' num2str(i_fish)]);
    
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
    
    %% [Next 2 cells] delete anatomically out-of-bound cells by hand
    % only execute once (manual choice 1):
    cHolder_Anat = []; % collect cIX of all out-of-bound cells
    
    %% (optional) step 1: set limit to only plot ventral layer, helps to find outliers within that layer
    % draw cells first
    % (choose start and stop index to draw; low numbers ~ ventral)
    I_start = 1;
    I_stop = 1*round(numcell_full/10);
    
    % plot those cells
    cIX = I_start:I_stop;
    gIX = round(cIX/1000)+1;
    figure('Position',[100 0 1300 900]);
    numK = round(cIX(end)/1000)+1;
    [~, dim_totimage] = BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);%,anat_zx,'hsv','full');
    
    % Here manually draw polygons around cells to discard.
    % Double click to connect last vertex to first vertex, then double click again within polygon to fix.
    % then click again to start drawing the next one
    % and when finished, break loop by pressing any key
    
    for i = 1:100, % NOMINAL LOOP, break manually
        h_poly_yx = impoly;
        wait(h_poly_yx); % double click to finalize position!
        % update finalized polygon in bright color
        setColor(h_poly_yx,[0 1 1]);
        
        A = sub2ind(size(anat_yx),CellXYZ(I_start:I_stop,1),CellXYZ(I_start:I_stop,2));
        MaskArray = createMask(h_poly_yx);
        B = find(MaskArray); % find indices of pixels within ROI
        cIX2 = find(ismember(A,B));
        cHolder_Anat = union(cHolder_Anat,cIX2);
        w = waitforbuttonpress;
        if w == 1, % press any key to break
            break; 
        end
    end
    
    %% step 2: plot all cells
    I_start = 1; % needs to be 1 for correct indexing
    I_stop = numcell_full;
    
    % ...here until the end of the cell is the exact dupliate of the last cell...
    cIX = I_start:I_stop;
    gIX = round(cIX/1000)+1;
    figure('Position',[100 0 1300 900]);
    numK = round(cIX(end)/1000)+1;
    [tot_image, dim_totimage] = BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);%,anat_zx,'hsv','full');
    
    % Here manually draw polygons around cells to discard.
    % Double click to connect last vertex to first vertex, then double click again within polygon to fix.
    % then click again to start drawing the next one
    % and when finished, break loop by pressing any key
    
    for i = 1:100, % NOMINAL LOOP, break manually
        h_poly_yx = impoly;
        wait(h_poly_yx); % double click to finalize position!
        % update finalized polygon in bright color
        setColor(h_poly_yx,[0 1 1]);
        
        A = sub2ind(size(anat_yx),CellXYZ(I_start:I_stop,1),CellXYZ(I_start:I_stop,2));
        MaskArray = createMask(h_poly_yx);
        B = find(MaskArray); % find indices of pixels within ROI
        cIX2 = find(ismember(A,B));
        cHolder_Anat = union(cHolder_Anat,cIX2);
        w = waitforbuttonpress;
        if w == 1, % press any key to break
            break; 
        end
    end
    
    %% find extra outliers on borders of image (can't always get with polygon)
    IX = find(CellXYZ(:,1)<15 | CellXYZ(:,1)> size(anat_stack,1)-15 ...
        | CellXYZ(:,2)<15 | CellXYZ(:,2)> size(anat_stack,2)-15);
    % for fish 9, eliminating bad plans: | CellXYZ(:,3)<4 | CellXYZ(:,3)> size(anat_stack,3)-3);
    cHolder_Anat = union(cHolder_Anat,IX);
    
    %% test Plot: all antomy outliers
    IX_inval_anat = cHolder_Anat; % rename
    cIX = IX_inval_anat;
    gIX = (1:length(cIX))';
    figure;
    BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);
    pause(2);
    
    %% test Plot: all remaining cells
    I_v_Holder = ones(1,numcell_full);
    I_v_Holder(IX_inval_anat) = 0;
    cIX = find(I_v_Holder);
    gIX = (1:length(cIX))';
    figure;
    BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);
    pause(2);
    
    %% ...and save
    temp = fullfile(data_dir,['Fish' num2str(i_fish) '_extrainfo_anat.mat']);
    save(temp,'IX_inval_anat');
    
    %%
    
    % save_masterdir = GetNestedDataDir();
    % save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    
    IX_inval = IX_inval_anat;
    filename = fullfile(save_dir,'OptionalInfo.mat');
    save(filename,'IX_inval','-append');
    
    filename = fullfile(save_dir,'AdditionalInfo.mat');
    save(filename,'IX_inval_anat','-append');
    disp('anat outliers saved');
    
end
