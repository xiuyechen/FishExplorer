function IX_m = GetMultifishClusIndex(IX,i_fish,range_fish)
currdir = GetCurrentDataDir();
load(fullfile(currdir,'BasicMultiFishInfo.mat','CellXYZ_multi'));
if ~ismember(i_fish,range_fish),
    disp('i_fish not included in fish range');
    IX_m = [];
    return;
end
ix = find(range_fish==i_fish);

numClus_multi = zeros(18,1);
for i = 1:18,
    numClus_multi(i) = max(gIX_Auto07_multi{i});
end


U = unique(gIX_Auto07_multi);
offset = length(cell2mat(CellXYZ_multi(range_fish(1:ix-1))));
IX_m = IX+offset;

end