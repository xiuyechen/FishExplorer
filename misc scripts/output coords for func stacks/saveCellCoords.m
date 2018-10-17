LeftRightFlag = gIX;
SelectCellXYZ_norm = CellXYZ_norm(cIX,:);
AllCellXYZ_norm = CellXYZ_norm;

save('CellCoords_norm.mat','AllCellXYZ_norm','LeftRightFlag','SelectCellXYZ_norm');

% save anat tiff stack