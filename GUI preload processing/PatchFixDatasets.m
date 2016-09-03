M_stimset = GetFishStimset();
M_dir = GetFishDirectories();
dataset_masterdir = GetNestedDataDir();

dim_norm = [1406,621,138];

%% 
range_fish = 12:18;%GetFishRange();

for i_fish = range_fish,
    disp(['i_fish = ', num2str(i_fish)]);

    %% load data
    dataset_dir = fullfile(dataset_masterdir,['subject_' num2str(i_fish)]);
%     s = load(fullfile(save_dir,'CoreInfo.mat'));
%     names = fieldnames(s);
    
    
    %%
    data_dir = M_dir{i_fish};
    if exist(fullfile(data_dir,'cell_resp_dim_lowcut.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim_lowcut.mat'));
    elseif exist(fullfile(data_dir,'cell_resp_dim.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim.mat'));
    else
        errordlg('find data to load!');
    end
    
    load(fullfile(data_dir,'cell_info.mat'));
        %% load anatomy
    tiffname = fullfile(data_dir,'ave.tif');
    info = imfinfo(tiffname,'tiff');
    nPlanes = length(info);
    s1 = info(1).Height;
    s2 = info(1).Width;
    ave_stack = zeros(s1,s2,nPlanes);
    for i=1:nPlanes,
        ave_stack(:,:,i) = imread(tiffname,i);
    end
    %%
     anat_stack = fliplr(ave_stack);
    [s1,s2,~] = size(anat_stack);
    for i_cell = 1:length(cell_info),
        % fix '.center'
        cell_info(i_cell).center(2) = s2-cell_info(i_cell).center(2)+1; %#ok<*SAGROW>
    end
    
    numcell_full = cell_resp_dim(1);
    temp = [cell_info(:).center];
    XY = reshape(temp',[],numcell_full)';
    Z = [cell_info.slice]';
    CellXYZ = horzcat(XY,Z);

    save(fullfile(dataset_dir,'CoreInfo.mat'),'CellXYZ','-append');
    
    %% flip!
    [CellXYZ_norm_raw,IX_inval_norm] = GetNormCellCord(i_fish);
    CellXYZ_norm = CellXYZ_norm_raw;
    CellXYZ_norm(:,2) = dim_norm(2) - CellXYZ_norm_raw(:,2) + 1;
    
    %%
    save_masterdir = GetCurrentDataDir();
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(save_dir,'data_full.mat'),'data');
    data.CellXYZ = CellXYZ;
    data.CellXYZ_norm = CellXYZ_norm;
%     names = fieldnames(data);
    
%     for i = 1:length(names), % use loop to save variables into fields of 'data'
%         eval(['data.',names{i},'=', names{i},';']);
%     end
    
    save(fullfile(save_dir,'data_full.mat'),'data');
end