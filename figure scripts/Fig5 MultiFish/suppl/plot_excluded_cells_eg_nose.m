% load 'IX_inval_anat' from 'Fish8_extrainfo_anat.mat' under 'Janelia2015'

absIX = getappdata(hfig,'absIX');
cIX_inval = IX_inval_anat(IX_inval_anat<6920);
cIX_val = setdiff(absIX,cIX_inval);
cIX_val_ds = cIX_val(1:1:end*0.1);
cIX = [cIX_val_ds;cIX_inval];
gIX = [ones(size(cIX_val_ds));2*ones(size(cIX_inval))];
setappdata(hfig,'clrmap_name','hsv_old');
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);

% HACK: go into function and change "cIX_abs = absIX(cIX);" to "cIX_abs = cIX;"
[h,im_full] = DrawCellsOnAnat(I);