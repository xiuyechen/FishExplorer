%%% Loading - Step 1. 
% Summary:
% load data from 'cell_resp.stackf' etc, and saved into
% 'FishX_direct_load.mat' (X = number) in the same directory.
% Stimulus parsing is done with function 'StimulusKeyPreprocessing', but
% one is advised to also inspect the results manually.

% Validation:
% There is an option 'isDiscard50' to discard the noisier 50%, disabled now.
% There is also a section to manually select (draw polygons) regions of
% cells to discard that are anatomical outliers. 

% Format:
% Cell-responses are converted into z-scores and kept that way through out.
% CInfo (cell info for the valid cells) stores the anatomical positions.
% For the rest of parameters saved, see the end of this file. 

% Next:
% For the following Loading steps, only this 'FishX_direct_load.mat' is 
% needed, the rest will be cleared from working memory.

%% 
clear all;close all;clc

code_dir = 'C:\Users\xiuye\Dropbox\Github\FishExplorer';
addpath(genpath(code_dir));
%% Set Manaully!

i_fish = 1; 

M_dir = GetFishDirectories();

M_tcutoff = {[],[],[],3000,3600,4000,1800,... % F1-7
    }

fpsec = 1.97; % Hz

%% load data
disp(['load data: fish ' num2str(i_fish)]);
datadir = M_dir{i_fish};
if exist(fullfile(datadir,'cell_resp_dim_lowcut.mat'), 'file'),
  load(fullfile(datadir,'cell_resp_dim_lowcut.mat'));
elseif exist(fullfile(datadir,'cell_resp_dim.mat'), 'file'),
    load(fullfile(datadir,'cell_resp_dim.mat'));
else
    errordlg('find data to load!');
end

load(fullfile(datadir,'cell_info.mat'));
if exist(fullfile(datadir,'cell_resp_lowcut.stackf'), 'file'),
  cell_resp_full = read_LSstack_fast_float(fullfile(datadir,'cell_resp_lowcut.stackf'),cell_resp_dim);
elseif exist(fullfile(datadir,'cell_resp.stackf'), 'file'),
    cell_resp_full = read_LSstack_fast_float(fullfile(datadir,'cell_resp.stackf'),cell_resp_dim);
else
    errordlg('find data to load!');
end

load(fullfile(datadir,'frame_turn.mat'));

%% load anatomy
tiffname = fullfile(datadir,'ave.tif');
info = imfinfo(tiffname,'tiff');
nPlanes = length(info);
s1 = info(1).Height;
s2 = info(1).Width;
ave_stack = zeros(s1,s2,nPlanes);
for i=1:nPlanes,
    ave_stack(:,:,i) = imread(fullfile(datadir,'ave.tif'),i);
end

% x-y view
im = max(ave_stack,[],3);
out=imNormalize99(im);
anat_yx = repmat(out,[1 1 3]);

% y-z view
im = squeeze(max(ave_stack,[],2));
out=imNormalize99(im);
anat_yz = repmat(out,[1 1 3]);

% x-z view
im = squeeze(max(ave_stack,[],1));
out=imNormalize99(im);
out = flipud(out'); %%%% empirically necessary...
anat_zx = repmat(out,[1 1 3]);

dimv_yx = size(anat_yx);
dimv_yz = size(anat_yz);
dimv_zx = size(anat_zx);

%% NEW: fix the left-right flip in the anatomy stack and subsequent cell_info
ave_stack = fliplr(ave_stack);
anat_yx = fliplr(anat_yx);
anat_zx = fliplr(anat_zx);

[s1,s2,~] = size(ave_stack);
tic
for i_cell = 1:length(cell_info),
    % fix '.center'
    cell_info(i_cell).center(2) = s2-cell_info(i_cell).center(2)+1; %#ok<*SAGROW>
    % fix '.inds'
    IX = cell_info(i_cell).inds;
    [I,J] = ind2sub([s1,s2],IX);
    J = s2-J+1;
    cell_info(i_cell).inds = sub2ind([s1,s2],I,J);
    % fix '.x_minmax'
    cell_info(i_cell).x_minmax(1) = s2-cell_info(i_cell).x_minmax(1)+1;
    cell_info(i_cell).x_minmax(2) = s2-cell_info(i_cell).x_minmax(2)+1;
end
toc

%% reformat coordinates
nCells = cell_resp_dim(1);
temp = [cell_info(:).center];
XY = reshape(temp',[],nCells)';
Z = [cell_info.slice]';
CellXYZ = horzcat(XY,Z);

%% Crop end of experiment when cell-segmentation drift > ~1um

%

%% Save mat files
% tic
% temp = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_full.mat']);
% varList = {'cell_resp_full','cell_resp_z_full','nCells','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
% % varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
% save(temp,varList{:},'-v7.3');
% toc

%% (if needed)
% save(temp,'somethingelse','-append');

poolobj=parpool(8);

range_fish = 1; %1:8

%
for i_fish = range_fish,
    disp(['i_fish = ', num2str(i_fish)]);
    
%     % loading
%     datadir = M_dir{i_fish};
%     file = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_nodiscard.mat']);
%     load(file,'CR_raw');
    %     varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn','periods'};
    
    %% detrend
    
    CR_dtr = zeros(size(cell_resp_full));
    CR_z_dtr = CR_dtr;
    tmax=size(cell_resp_full,2);
    nCells=size(cell_resp_full,1);
    
    tic
    parfor i=1:nCells,
        cr = cell_resp_full(i,:);
        crd = 0*cr;
        for j=1:100:tmax,
            if j<=150,
                tlim1 = 1;
                tlim2 = 300;
            elseif j>tmax-150,
                tlim1 = tmax-300;
                tlim2 = tmax;
            else
                tlim1 = j-150;
                tlim2 = j+150;
            end
            crr = cr(tlim1:tlim2);
            crd(max(1,j-50):min(tmax,j+50)) = prctile(crr,15);
        end
        if mod(i,100)==0,
            disp(num2str(i));
        end
        CR_dtr(i,:) = cr-crd;
        CR_z_dtr(i,:) = zscore(cr-crd);
    end
    toc
    
    CR_dtr = single(CR_dtr);
    CR_z_dtr = single(CR_z_dtr);
    
    %% save into .mat
    %     save(file,'CR_dtr','-append');
    %% Save mat files
    tic
    temp = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_detrend.mat']);
    varList = {'CR_dtr','CR_z_dtr','nCells','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
    % varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
    save(temp,varList{:},'-v7.3');
    toc

end

%%
delete(poolobj);

%% next run step 2: detrend


