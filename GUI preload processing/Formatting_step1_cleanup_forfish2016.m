% Formatting_step1_cleanup

clear all;close all;clc

code_dir = 'C:\Users\Xiu\Dropbox\FishExplorer';
addpath(genpath(code_dir));

%% Set Manaully!
% file directories
M_dir = GetFishDirectories();
M_tcutoff = {3500,[],[],3500,3600,4000,1800,[],... % Fish 1-8
    [],5000,[],7000,6500,4500, ... % fish 9-14 % fish 12 last ran with []! 1/22/16
    [],2200,[],7000}; % fish 15-18 with GAD in red channel

M_fpsec = {1.97,1.97,1.97,1.97,1.97,1.97,1.97,...% in Hz % Fish 1-7
    1.97,1.93,1.97,1.84,... % Fish 8-11 % Actually Fish 10 freq unrecorded
    2.56,2.33,2.27,... % Fish 12-14
    2.33,2.33,2.38,2.38}; % Fish 15-18

% also check manual frame-correction section below

%%
% poolobj=parpool(8);
%%
range_fish = GetFishRange();

for i_fish = range_fish,
    tic
    disp(['i_fish = ', num2str(i_fish)]);    

    fpsec = M_fpsec{i_fish};
    
    %% load data
    disp(['load data: fish ' num2str(i_fish)]);
    data_dir = M_dir{i_fish};
    if exist(fullfile(data_dir,'cell_resp_dim_lowcut.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim_lowcut.mat'));
    elseif exist(fullfile(data_dir,'cell_resp_dim.mat'), 'file'),
        load(fullfile(data_dir,'cell_resp_dim.mat'));
    else
        errordlg('find data to load!');
    end
    
    load(fullfile(data_dir,'cell_info.mat'));

    if exist(fullfile(data_dir,'frame_turn_new.mat'), 'file'),
        load(fullfile(data_dir,'frame_turn_new.mat'));
    elseif exist(fullfile(data_dir,'frame_turn.mat'), 'file'),
        load(fullfile(data_dir,'frame_turn.mat'));
    end
    frame_keys = frame_turn;
    
    if i_fish<12,
        if exist(fullfile(data_dir,'cell_resp_lowcut.stackf'), 'file'),
            TimeSeries = read_LSstack_fast_float(fullfile(data_dir,'cell_resp_lowcut.stackf'),cell_resp_dim);
        elseif exist(fullfile(data_dir,'cell_resp.stackf'), 'file'),
            TimeSeries = read_LSstack_fast_float(fullfile(data_dir,'cell_resp.stackf'),cell_resp_dim);
        else
            errordlg('find data to load!');
        end
    else
        TimeSeries = read_LSstack_fast_float(fullfile(data_dir,'CR_detrend.stackf'),cell_resp_dim);
    end
    
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
    
    %% fix the left-right flip in the anatomy stack and subsequent cell_info
    anat_stack = fliplr(ave_stack);
    %     anat_stack(:,113:1092,:) = anat_stackCopy;
    
    [s1,s2,~] = size(anat_stack);
    for i_cell = 1:length(cell_info),
        % fix '.center'
        cell_info(i_cell).center(2) = s2-cell_info(i_cell).center(2)+1; %#ok<*SAGROW>
    end

    %% reformat coordinates
    numcell_full = cell_resp_dim(1);
    temp = [cell_info(:).center];
    XY = reshape(temp',[],numcell_full)';
    Z = [cell_info.slice]';
    CellXYZ = horzcat(XY,Z);
    
        
%     %% Manual occasional frame corrections: ONLY RUN ONCE!!!!!!!!!
%         
%     if true,
%         % add 1 frame at start: ------------- for which fish again? all till F11?
%         cell_resp_full = horzcat(cell_resp_full(:,1),cell_resp_full);
%         frame_keys = vertcat(frame_keys(1,:),frame_keys);
%     end
% 
%     % correction of an error in Fish #1
%     if i_fish == 1,
%         cell_resp_full = horzcat(cell_resp_full(:,1:20),cell_resp_full);
%     end
%     
    %% Crop end of experiment (when cell-segmentation drift > ~1um, manually determined)
    if i_fish==13,
        TimeSeries = TimeSeries(:,305:end);
    end
    
    if ~isempty(M_tcutoff{i_fish}),
        TimeSeries = TimeSeries(:,1:M_tcutoff{i_fish});
%         cell_resp_full = cell_resp_full(:,1:M_tcutoff{i_fish});
        frame_keys = frame_keys(1:M_tcutoff{i_fish},:);
    end

    %% detrend    
%     CR_dtr = zeros(size(cell_resp_full));
%     tmax=size(cell_resp_full,2);
%     numcell_full=size(cell_resp_full,1);
% 
%     parfor i=1:numcell_full,
%         cr = cell_resp_full(i,:);
%         crd = 0*cr;
%         for j=1:100:tmax,
%             if j<=150,
%                 tlim1 = 1;
%                 tlim2 = 300;
%             elseif j>tmax-150,
%                 tlim1 = tmax-300;
%                 tlim2 = tmax;
%             else
%                 tlim1 = j-150;
%                 tlim2 = j+150;
%             end
%             crr = cr(tlim1:tlim2);
%             crd(max(1,j-50):min(tmax,j+50)) = prctile(crr,15);
%         end
% %         if mod(i,100)==0,
% %             disp(num2str(i));
% %         end
%         CR_dtr(i,:) = cr-crd;
%     end
%         
%     TimeSeries = single(CR_dtr);

    
    %% registered coordinates from morphing to ZBrain
    if i_fish<=18,
        [CellXYZ_norm,IX_inval_norm] = GetNormCellCord(i_fish);
    else % place holder
        CellXYZ_norm = CellXYZ;
        IX_inval_norm = [];
    end    
    
    %% Save files
tic
    % directory to save data formatted for distribution:
    save_masterdir = GetNestedDataDir();
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % save time-series
    h5create(fullfile(save_dir,'TimeSeries.h5'),'/TimeSeries',size(TimeSeries),'Datatype','single','ChunkSize',[1000 100]);
    h5write(fullfile(save_dir,'TimeSeries.h5'),'/TimeSeries',TimeSeries);

    % save .mat files
    filename = fullfile(save_dir,'CoreInfo.mat');
    varList = {'CellXYZ','anat_stack','fpsec'};
    save(filename,varList{:});
    
    filename = fullfile(save_dir,'OptionalInfo.mat');
    varList = {'numcell_full','CellXYZ_norm','IX_inval_norm'};
    save(filename,varList{:});    
    
    filename = fullfile(save_dir,'AdditionalInfo.mat');
    save(filename,'frame_keys');
toc
    %%
    toc
end

%%
% delete(poolobj);

