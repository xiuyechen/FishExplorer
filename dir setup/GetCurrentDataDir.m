function data_dir = GetCurrentDataDir()
if ispc,
    data_dir = 'C:\Janelia2015\GUI_data';
elseif isunix
    data_dir = '/home/xiuye/Research/GUI_data';
end