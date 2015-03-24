%% 
clear all;close all;clc

i_fish = 9; %%%% <------- manual set

%% constants
M_dir = {'F:\Janelia2014\Fish1_16states_30frames';
    'F:\Janelia2014\Fish2_20140714_2_4_16states_10frames';
    'F:\Janelia2014\Fish3_20140715_1_1_16_states_10frames';
    'F:\Janelia2014\Fish4_20140721_1_8_16states_20frames';
    'F:\Janelia2014\Fish5_20140722_1_2_16states_30frames';
    'F:\Janelia2014\Fish6_20140722_1_1_3states_30,40frames';
    'F:\Janelia2014\Fish7_20140722_2_3_3states_30,50frames';
    'F:\Janelia2014\Fish8_20141222_2_2_7d_PT_3OMR_shock_lowcut';
    'F:\Janelia2014\Fish9_20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356'};
    
global fpsec;
fpsec = 1.97; % Hz

M_period = {480,160,160,320,480,280,300,{120,150,360},{1000}};
periods = M_period{i_fish};

%% load data
datadir = M_dir{i_fish};
% load(fullfile(datadir,'cell_resp_dim_lowcut.mat'));
load(fullfile(datadir,'cell_resp_dim.mat'));
load(fullfile(datadir,'cell_info.mat'));
cell_resp = read_LSstack_fast_float(fullfile(datadir,'cell_resp.stackf'),cell_resp_dim);
load(fullfile(datadir,'frame_turn.mat'));

% frame_turn adjustments
frame_turn = frame_turn';
if i_fish==6 || i_fish==7,
    temp = round(frame_turn(17,:));
    temp(temp==0)=-3;temp(temp==2)=-3;temp(temp==3)=-2; temp(temp==-3)=3;temp(temp==-2)=2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    frame_turn(17,:)=temp;
end

% shift/photostate etc is only prelim here, also rep avr is prelim
% (will be included in step2 and saved into CONST, not direct_load)
M_period_0 = {480,320,320,320,480,280,300,400,500};
period_0 = M_period_0{i_fish};
if i_fish == 1,
    shift = 460;
else
    shift = period_0;
end
photostate = round(frame_turn(17,1+shift:period_0+shift));

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
% anatomy_image_yx = repmat(imNormalize99(max(ave_stack,[],3)),[1 1 3]);
% anatomy_image_yz = repmat(imNormalize99(squeeze(max(ave_stack,[],2))),[1 1 3]);

% imNormalize99.m
% x-y view
im = max(ave_stack,[],3);
im=double(im);
temp=sort(im(:),'descend');
th1=temp(round(length(im(:))/100));
th2=min(im(:));
out=(im-th2)/(th1-th2);
out(out>1)=1;
anat_yx = repmat(out,[1 1 3]);

% y-z view
im = squeeze(max(ave_stack,[],2));
im=double(im);
temp=sort(im(:),'descend');
th1=temp(round(length(im(:))/100));
th2=min(im(:));
out=(im-th2)/(th1-th2);
out(out>1)=1;
anat_yz = repmat(out,[1 1 3]);

% x-z view
im = squeeze(max(ave_stack,[],1));
im=double(im);
temp=sort(im(:),'descend');
th1=temp(round(length(im(:))/100));
th2=min(im(:));
out=(im-th2)/(th1-th2);
out(out>1)=1;
out = flipud(out');
anat_zx = repmat(out,[1 1 3]);

dimv_yx = size(anat_yx);
dimv_yz = size(anat_yz);
dimv_zx = size(anat_zx);

%% crop, average
numcell=length(cell_info);
totlen=cell_resp_dim(2);
nrep=floor(totlen/period_0)-1;

skiplist=[];
IX_rep=setdiff(1:nrep, skiplist);
IX=zeros(period_0*length(IX_rep),1);

for i=1:length(IX_rep)
    IX(period_0*(i-1)+1:period_0*i)=period_0*(IX_rep(i)-1)+1+shift:period_0*(IX_rep(i))+shift;
end

% cell_resp=cell_resp(:,1:nrep*period);
cell_resp_ave=mean(reshape(cell_resp(:,IX),[numcell period_0 length(IX_rep)]),3);
cell_resp_ave_z = zscore(cell_resp_ave')';
cell_resp_z = zscore(cell_resp')';

%%
% tiffName = fullfile(datadir,['Fish' num2str(i_fish) '_allcells.tif']);
% WriteZstack_LSh(1:numcell,ones(1,numcell),ave_stack,cell_info,tiffName) % (cIX,gIX,ave_stack,cell_info,tiffName,cIX_0)

%% Validating all cells

%%%% Round 1: discard 50% noisy cells based on std of zscore of baseline 

% Compute std of dimest 10% of frames for each cell. 
prc = prctile(cell_resp_z,10,2);
STD_full = zeros(numcell,1);
for i = 1:numcell,
    ix = find(cell_resp_z(i,:)<prc(i));
    STD_full(i) = std(cell_resp_z(i,ix));
end
% Set threshold at 50%, i.e. discard 50% of all cells
thr = prctile(STD_full,50); 

temp = find(STD_full>thr);
[~,I] = sort(STD_full(temp));
cIX = temp(I);

% visualize cells to discard
% gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
gIX = round(cIX/1000)+1; % option: view anatomical z index
M = cell_resp_ave_z(cIX,:);
BasicPlotMaps(cIX,gIX,M,cell_info,photostate,anat_yx,anat_yz,anat_zx);

% I_v: index of valid cells
I_v_Holder = ones(1,numcell); 
I_v_Holder(cIX) = 0;
I_v = find(I_v_Holder);

CRAZ = cell_resp_ave_z(I_v,:);
STD = STD_full(I_v);
nCells = length(I_v);
CInfo = cell_info(I_v);

%% visualize valid cells, sorted by noise
[~,I] = sort(STD);
cIX = I;
gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
M = CRAZ(I,:);
BasicPlotMaps(cIX,gIX,M,CInfo,photostate,anat_yx,anat_yz,anat_zx);

%%
%%%% Round 2: delete anatomically out-of-bound cells by hand

% draw cells first
% manual choice 1:
I1 = 1;
I2 = 2000; % manually set limit to see the layers better

%% manual choice 2:
I1 = 1;
I2 = length(I_v);
%%
cIX = I1:I2;
gIX = round(cIX/1000)+1;
figure('Position',[100 0 1300 900]);
h2=axes('Position',[0.61, 0.06, 0.37, 0.85]); % ~subplot
numK = round(cIX(end)/1000)+1;
[tot_image, dim_totimage] = DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

%% hand draw ROI group masks
% only execute once (manual choice 1):
cHolder_Anat = []; % collect cIX of all out-of-bound cells
%% (break with Ctrl+C)
MaskArray = zeros(dim_totimage(1), dim_totimage(2));
k_zres = 20;
% draw polygon around outliers on Y-X view
for i = 1:10, % NOMINAL LOOP, break by pressing key instead of button
    h_poly_yx = impoly;
    wait(h_poly_yx); % double click to finalize position!
    % update finalized polygon in bright color
    setColor(h_poly_yx,[0 1 1]);

    IJs = reshape([CInfo.center],2,[])';    
    A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
    MaskArray = createMask(h_poly_yx);
    MaskArray(1:dimv_zx*k_zres+10,:) = [];
    B = find(MaskArray); % find indices of pixels within ROI
    cIX2 = find(ismember(A,B));
    cHolder_Anat = union(cHolder_Anat,cIX2);
    w = waitforbuttonpress;
    if w == 1,
        break;
    end
end
% end of manual choices loop

%% find extra's on top/bottom border of image (can't get with polygon)
temp = [CInfo(:).center];
XY = reshape(temp',[],nCells)';
IX = find(XY(:,1)<8 | XY(:,1)> 2040);

cHolder_Anat = union(cHolder_Anat,IX);

%% test Plot traces + map for antomy outliers
cIX = cHolder_Anat;
gIX = (1:length(cIX))';
M = CRAZ(cIX,:);
% BasicPlotMaps(cIX,gIX,M,CInfo,photostate,anat_yx,anat_yz,anat_zx);
figure;
DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

%% Clean up anatomy outliers
cIX_Invalid_Anat = I_v(cHolder_Anat); % convert back to index for original full set
I_v_Holder(cIX_Invalid_Anat) = 0;
I_v2 = find(I_v_Holder);
%%
CR_raw = cell_resp(I_v2,:);
% STD_10prc = STD_full(I_v);
nCells = length(I_v2);
CInfo = cell_info(I_v2);

%% Plot all valid cells to save
cIX = (1:nCells)';
gIX = round(cIX/1000)+1;
CRAZ = cell_resp_ave_z(I_v2,:);
M = CRAZ(cIX,:);
figure;
DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,'hsv','full');

%% Save mat files
CInfo_full = cell_info; % save in next round! not saved in current mats!!!!!!

temp = fullfile(datadir,['Fish' num2str(i_fish) '_full_extrainfo.mat']);
save(temp,'cHolder_Anat','cIX_Invalid_Anat','I_v','STD_full','CInfo_full'); % CInfo_full not saved in current mats!!

temp = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load.mat']);
varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn','periods'};
save(temp,varList{:});

%%
% save(temp,'anat_zx','-append');

%% next run Detrend_all, before running step 2


