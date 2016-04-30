M_dir = GetFishDirectories();

Coords7fish = cell(1,7);
range_fish = 12:18,
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    data_dir = M_dir{i_fish};
    if exist(fullfile(data_dir,'cell_resp_dim_lowcut.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim_lowcut.mat'));
    elseif exist(fullfile(data_dir,'cell_resp_dim.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim.mat'));
    end
    load(fullfile(data_dir,'cell_info.mat'));
    
        numcell_full = cell_resp_dim(1);
    temp = [cell_info(:).center];
    XY = reshape(temp',[],numcell_full)';
    Z = [cell_info.slice]';
    CellXYZ = horzcat(XY,Z);
    Coords7fish{i} = CellXYZ;
%     LoadFishDataWithoutTS(hfig,i_fish);
%     Coords7fish{i_fish} = getappdata(hfig,'CellXYZ');
end

data_dir = GetCurrentDataDir();
save(fullfile(data_dir,'Coords7fish.mat'),'Coords7fish');