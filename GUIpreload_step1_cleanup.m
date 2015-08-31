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

%% Set Manaully!

i_fish = 13; 

M_dir = GetFishDirectories();

fpsec = 1.97; % Hz

%% load data
disp(['load data: fish ' num2str(i_fish)]);
datadir = M_dir{i_fish};
if exist(fullfile(datadir,'cell_resp_dim.mat'), 'file'),
  load(fullfile(datadir,'cell_resp_dim.mat'));
elseif exist(fullfile(datadir,'cell_resp_dim_lowcut.mat'), 'file'),
    load(fullfile(datadir,'cell_resp_dim_lowcut.mat'));
else
    errordlg('find data to load!');
end

load(fullfile(datadir,'cell_info.mat'));
if exist(fullfile(datadir,'cell_resp.stackf'), 'file'),
  cell_resp_full = read_LSstack_fast_float(fullfile(datadir,'cell_resp.stackf'),cell_resp_dim);
elseif exist(fullfile(datadir,'cell_resp_lowcut.stackf'), 'file'),
    cell_resp_full = read_LSstack_fast_float(fullfile(datadir,'cell_resp_lowcut.stackf'),cell_resp_dim);
else
    errordlg('find data to load!');
end

load(fullfile(datadir,'frame_turn.mat'));

% parse stimulus - should check manually to be sure!
[stimset,stim] = StimulusKeyPreprocessing(frame_turn,i_fish);

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
    cell_info(i_cell).center(2) = s2-cell_info(i_cell).center(2)+1;
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

%% compute z-score
disp('compute z-score')
cell_resp_z_full = zscore(cell_resp_full')';
nCells = cell_resp_dim(1);

%% Validating all cells

isDiscard50 = false; % option to pre-screen
if isDiscard50,
    %% Round 1: discard 50% noisy cells based on std of zscore of baseline
    
    % Compute std of dimest 10% of frames for each cell.
    prc = prctile(cell_resp_z_full,10,2);
    STD_full = zeros(nCells,1);
    for i = 1:nCells,
        ix = find(cell_resp_z_full(i,:)<prc(i));
        STD_full(i) = std(cell_resp_z_full(i,ix));
    end
    % Set threshold at 50%, i.e. discard 50% of all cells
    thr = prctile(STD_full,50);
    
    temp = find(STD_full>thr);
    [~,I] = sort(STD_full(temp));
    cIX = temp(I);
    
    % visualize cells to discard
    % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
    gIX = round(cIX/1000)+1; % option: view anatomical z index
    M = cell_resp_z_full(cIX,1:1000);
    BasicPlotMaps(cIX,gIX,M,cell_info,stim,anat_yx,anat_yz,anat_zx);
    
    % I_v: index of valid cells
    I_v_Holder = ones(1,nCells);
    I_v_Holder(cIX) = 0;
    I_v = find(I_v_Holder);
    
    cell_resp_z = cell_resp_z_full(I_v,:);
    STD = STD_full(I_v);
    nCells = length(I_v);
    CInfo_0 = cell_info(I_v);
    
    %% visualize valid cells, sorted by noise
    [~,I] = sort(STD);
    cIX = I;
    gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
    M = cell_resp_z(I,1:1000);
    BasicPlotMaps(cIX,gIX,M,CInfo_0,stim,anat_yx,anat_yz,anat_zx);

else    
    %% Round 1 alternative: full set
    cell_resp_z = cell_resp_z_full;
    I_v_Holder = ones(1,nCells);
    I_v = (1:nCells)'; % valid Index. here not actually screened, I_v is fullset

    cIX = I_v; % 'cell-Index'
    gIX = (round((1:length(cIX))/1000)+1)'; % sorted 'group-Index'
    M = cell_resp_z_full(cIX,1:1000);
    CInfo_0 = cell_info;
    BasicPlotMaps(cIX,gIX,M,CInfo_0,stim,anat_yx,anat_yz,anat_zx);
end

%% %%%% Round 2: delete anatomically out-of-bound cells by hand
% only execute once (manual choice 1):
cHolder_Anat = []; % collect cIX of all out-of-bound cells


% draw cells first
% (choose start and stop index to draw; low numbers ~ ventral)

%% (optional) step 1: set limit to only plot ventral layer, helps to find outliers within that layer
I_start = 1;
I_stop = round(length(I_v)/10);

% plot those cells
cIX = I_start:I_stop;
gIX = round(cIX/1000)+1;
figure('Position',[100 0 1300 900]);
numK = round(cIX(end)/1000)+1;
[tot_image, dim_totimage] = DrawClustersOnMap_LSh(CInfo_0,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

% Here manually draw polygons around cells to discard. 
% Double click to connect last vertex to first vertex, then double click again within polygon to fix.
% then click again to start drawing the next one
% and when finished, break loop with Ctrl+C...
MaskArray = zeros(dim_totimage(1), dim_totimage(2));
k_zres = 20;
% draw polygon around outliers on Y-X view
for i = 1:100, % NOMINAL LOOP, break manually (with Ctrl+C, somehow the design didn't work)
    h_poly_yx = impoly;
    wait(h_poly_yx); % double click to finalize position!
    % update finalized polygon in bright color
    setColor(h_poly_yx,[0 1 1]);

    IJs = reshape([CInfo_0(I_start:I_stop).center],2,[])';    
    A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
    MaskArray = createMask(h_poly_yx);
    MaskArray(1:dimv_zx*k_zres+10,:) = [];
    B = find(MaskArray); % find indices of pixels within ROI
    cIX2 = find(ismember(A,B));
    cHolder_Anat = union(cHolder_Anat,cIX2);
    w = waitforbuttonpress;
    if w == 1,
        break; % this doesn't work ?!
    end
end

%% step 2: plot all cells (including the ones in the ventral layer, don't 
I_start = 1;
I_stop = length(I_v); 

% ...here until the end of the cell is the exact dupliate of the last cell...
cIX = I_start:I_stop;
gIX = round(cIX/1000)+1;
figure('Position',[100 0 1300 900]);
numK = round(cIX(end)/1000)+1;
[tot_image, dim_totimage] = DrawClustersOnMap_LSh(CInfo_0,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

% Here manually draw polygons around cells to discard. 
% Double click to connect last vertex to first vertex, then double click again within polygon to fix.
% then click again to start drawing the next one
% and when finished, break loop with Ctrl+C...
MaskArray = zeros(dim_totimage(1), dim_totimage(2));
k_zres = 20;
% draw polygon around outliers on Y-X view
for i = 1:100, % NOMINAL LOOP, break manually (with Ctrl+C, somehow the design didn't work)
    h_poly_yx = impoly;
    wait(h_poly_yx); % double click to finalize position!
    % update finalized polygon in bright color
    setColor(h_poly_yx,[0 1 1]);

    IJs = reshape([CInfo_0(I_start:I_stop).center],2,[])';    
    A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
    MaskArray = createMask(h_poly_yx);
    MaskArray(1:dimv_zx*k_zres+10,:) = [];
    B = find(MaskArray); % find indices of pixels within ROI
    cIX2 = find(ismember(A,B));
    cHolder_Anat = union(cHolder_Anat,cIX2);
    w = waitforbuttonpress;
    if w == 1,
        break; % this doesn't work ?!
    end
end

%% find extra outliers on top/bottom border of image (can't always get with polygon)
temp = [CInfo_0(:).center];
XY = reshape(temp',[],nCells)';
IX = find(XY(:,1)<8 | XY(:,1)> 2040);

cHolder_Anat = union(cHolder_Anat,IX);

%% test Plot: all antomy outliers
cIX = cHolder_Anat;
gIX = (1:length(cIX))';
figure;
DrawClustersOnMap_LSh(CInfo_0,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

%% Clean up anatomy outliers
cIX_Invalid_Anat = I_v(cHolder_Anat); % convert back to index for original full set
I_v_Holder(cIX_Invalid_Anat) = 0;
I_v2 = find(I_v_Holder);

CR_raw = cell_resp_z_full(I_v2,:);
nCells = length(I_v2);
CInfo = cell_info(I_v2);

%% Plot all valid cells to save
cIX = (1:nCells)';
gIX = round(cIX/1000)+1;
M = CR_raw(cIX,1:1000);
figure;
DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

%% Save mat files
CInfo_full = cell_info; 
tic
temp = fullfile(datadir,['Fish' num2str(i_fish) '_full_extrainfo.mat']);
save(temp,'cHolder_Anat','cIX_Invalid_Anat','I_v','I_v2','CInfo_full','-v7.3'); % CInfo_full not saved in mats before 4/24/15!! %% 'STD_full' saved before

temp = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_nodiscard.mat']);
varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
save(temp,varList{:},'-v7.3');
toc

%% (if needed)
% save(temp,'somethingelse','-append');

%% next run step 2: detrend


