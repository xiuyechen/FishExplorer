function IX_m = GetMultifishCellIndex(IX,i_fish)
currdir = GetCurrentDataDir();
load(fullfile(currdir,'BasicMultiFishInfo.mat','CellXYZ_multi'));

offset = length(cell2mat(CellXYZ_multi(1:i_fish-1)));
IX_m = IX+offset;

end