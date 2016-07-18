function LoadFullFish(hfig,i_fish,isFullData,hdf5_dir,mat_dir)
M_fish_set = getappdata(hfig,'M_fish_set');
fishset = M_fish_set(i_fish);
setappdata(hfig,'fishset',fishset);

%% load data from file
disp(['loading fish ' num2str(i_fish) '...']);
loadStart = tic;
isLoadTS = 1;
if ~exist('hdf5_dir','var') || ~exist('mat_dir','var'),
    data_masterdir = getappdata(hfig,'data_masterdir');
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    if exist('isFullData','var'),
        if isFullData == -1,
            hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
            absIX = h5read(hdf5_dir,'/absIX');
            setappdata(hfig,'absIX',absIX);
            isLoadTS = 0; % don't load Time series                       
        elseif isFullData, % true
             hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
        else % false
            hdf5_dir = fullfile(data_dir,'TimeSeries_half.h5');
        end
    else % default without isfulldata input is full-data
        hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
    end
    mat_dir = fullfile(data_dir,'data_full.mat');
end

% load 'data'
load(mat_dir,'data'); % struct with many fields
names = fieldnames(data); % cell of strings
for i = 1:length(names),
    setappdata(hfig,names{i},eval(['data.',names{i}]));
end

if isLoadTS,
    % load time series (hdf5 file)
    CellResp = h5read(hdf5_dir,'/CellResp');
    CellRespZ = h5read(hdf5_dir,'/CellRespZ');
    CellRespAvr = h5read(hdf5_dir,'/CellRespAvr');
    CellRespAvrZ = h5read(hdf5_dir,'/CellRespAvrZ');
    absIX = h5read(hdf5_dir,'/absIX');
    
    setappdata(hfig,'CellResp',CellResp);
    setappdata(hfig,'CellRespZ',CellRespZ);
    setappdata(hfig,'CellRespAvr',CellRespAvr);
    setappdata(hfig,'CellRespAvrZ',CellRespAvrZ);
    setappdata(hfig,'absIX',absIX);
end

loadTime = toc(loadStart);
disp([num2str(loadTime) ' sec']);

end
