function LoadFishDataWithoutTS(hfig,i_fish,hdf5_dir,mat_dir)
disp(['loading fish #' num2str(i_fish) '...']);

M_fish_set = getappdata(hfig,'M_fish_set');
fishset = M_fish_set(i_fish);
setappdata(hfig,'fishset',fishset);

%% load data from file

disp(['load fish ' num2str(i_fish) '...']);

if ~exist('hdf5_dir','var') || ~exist('mat_dir','var'),
    data_masterdir = getappdata(hfig,'data_masterdir');
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
    mat_dir = fullfile(data_dir,'data_full.mat');
end

tic

% load 'data'
load(mat_dir,'data'); % struct with many fields
names = fieldnames(data); % cell of strings
for i = 1:length(names),
    setappdata(hfig,names{i},eval(['data.',names{i}]));
end

% load time series (hdf5 file)
absIX = h5read(hdf5_dir,'/absIX');

setappdata(hfig,'absIX',absIX);

toc

% setappdata(hfig,'isfullfish',1);
end
