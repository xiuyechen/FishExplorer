clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_dir = GetCurrentDataDir();
% save_dir = GetCurrentDataDir();
i_fish = 1;

%% load data
disp(['load fish ' num2str(i_fish)]);

M_dir = GetFishDirectories();
varList = {'numcell_full','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
load(fullfile(M_dir{i_fish},['Fish' num2str(i_fish) '_direct_load_full_v6.mat']),varList{:}); 

%% [Next 2 cells] delete anatomically out-of-bound cells by hand
% only execute once (manual choice 1):
cHolder_Anat = []; % collect cIX of all out-of-bound cells

%% (optional) step 1: set limit to only plot ventral layer, helps to find outliers within that layer
% draw cells first
% (choose start and stop index to draw; low numbers ~ ventral)
I_start = 1;
I_stop = round(numcell_full/10);

% plot those cells
cIX = I_start:I_stop;
gIX = round(cIX/1000)+1;
figure('Position',[100 0 1300 900]);
numK = round(cIX(end)/1000)+1;
[~, dim_totimage] = BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);%,anat_zx,'hsv','full');

% Here manually draw polygons around cells to discard. 
% Double click to connect last vertex to first vertex, then double click again within polygon to fix.
% then click again to start drawing the next one
% and when finished, break loop with Ctrl+C...

% MaskArray = zeros(dim_totimage(1), dim_totimage(2));
% k_zres = 20;
% draw polygon around outliers on Y-X view
for i = 1:100, % NOMINAL LOOP, break manually (with Ctrl+C, somehow the design didn't work)
    h_poly_yx = impoly;
    wait(h_poly_yx); % double click to finalize position!
    % update finalized polygon in bright color
    setColor(h_poly_yx,[0 1 1]);

    IJs = reshape([cell_info(I_start:I_stop).center],2,[])';    
    A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
    MaskArray = createMask(h_poly_yx);
%     MaskArray(1:dimv_zx*k_zres+10,:) = [];
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
I_stop = numcell_full; 

% ...here until the end of the cell is the exact dupliate of the last cell...
cIX = I_start:I_stop;
gIX = round(cIX/1000)+1;
figure('Position',[100 0 1300 900]);
numK = round(cIX(end)/1000)+1;
[tot_image, dim_totimage] = BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);%,anat_zx,'hsv','full');

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

    IJs = reshape([cell_info(I_start:I_stop).center],2,[])';    
    A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
    MaskArray = createMask(h_poly_yx);
%     MaskArray(1:dimv_zx*k_zres+10,:) = [];
    B = find(MaskArray); % find indices of pixels within ROI
    cIX2 = find(ismember(A,B));
    cHolder_Anat = union(cHolder_Anat,cIX2);
    w = waitforbuttonpress;
    if w == 1,
        break; % this doesn't work ?!
    end
end

%% find extra outliers on borders of image (can't always get with polygon)
% temp = [cell_info(:).center];
% XY = reshape(temp',[],numcell_full)';
IX = find(XY(:,1)<15 | XY(:,1)> s1-15 ...
        | XY(:,2)<15 | XY(:,2)> s2-15);
cHolder_Anat = union(cHolder_Anat,IX);

%% test Plot: all antomy outliers
IX_inval_anat = cHolder_Anat; % rename
cIX = IX_inval_anat;
gIX = (1:length(cIX))';
figure;
BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);

%% ...and save
temp = fullfile(datadir,['Fish' num2str(i_fish) '_extrainfo_anat.mat']);
save(temp,'IX_inval_anat');

%% test Plot: all remaining cells
I_v_Holder = ones(1,numcell_full);
I_v_Holder(IX_inval_anat) = 0;
cIX = find(I_v_Holder);
gIX = (1:length(cIX))';
figure;
BasicDrawCellsOnAnatProj(CellXYZ,cIX,gIX,numK,anat_yx,anat_yz);