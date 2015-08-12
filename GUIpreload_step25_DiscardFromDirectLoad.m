
clear all;close all;clc

M_dir = {'E:\Janelia2014\Fish1_16states_30frames';
    'E:\Janelia2014\Fish2_20140714_2_4_16states_10frames';
    'E:\Janelia2014\Fish3_20140715_1_1_16_states_10frames';
    'E:\Janelia2014\Fish4_20140721_1_8_16states_20frames';
    'E:\Janelia2014\Fish5_20140722_1_2_16states_30frames';
    'E:\Janelia2014\Fish6_20140722_1_1_3states_30,40frames';
    'E:\Janelia2014\Fish7_20140722_2_3_3states_30,50frames';
    'E:\Janelia2014\Fish8_20141222_2_2_7d_PT_3OMR_shock_lowcut';
    'E:\Janelia2014\Fish9_20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356';
    'E:\Janelia2014\Fish10_20150120_2_2_photo_OMR_prey_blob_blue_cy74_6d_20150120_231917'};

varList = {'CR_dtr','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};

i_fish = 10;

%% load data
disp(['load fish ' num2str(i_fish)]);
tic
datadir = M_dir{i_fish};
load(fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_nodiscard.mat']),varList{:});
toc

CR_dtr_full = CR_dtr;
CInfo_full = CInfo;
%% Test: discard 50% noisy cells based on std of zscore of baseline

% Compute std of dimest 10% of frames for each cell.
tic
prc = prctile(CR_dtr_full,10,2);
nCells = size(CR_dtr_full,1);
STD_full = zeros(nCells,1);
for i = 1:nCells,
    ix = find(CR_dtr_full(i,:)<prc(i));
    STD_full(i) = std(CR_dtr_full(i,ix));
end
toc
beep

%% Set threshold at 50%, i.e. discard 50% of all cells
perc_keep = 75;
thr = prctile(STD_full,perc_keep); 

temp = find(STD_full>thr);
[~,I] = sort(STD_full(temp));
cIX = temp(I);

%% visualize cells to discard
% gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
gIX = round(cIX/1000)+1; % option: view anatomical z index
M = CR_dtr_full(cIX,1:1000);
[~,stim] = StimulusKeyPreprocessing(frame_turn,i_fish);
BasicPlotMaps(cIX,gIX,M,CInfo,stim,anat_yx,anat_yz,anat_zx);

% I_v: index of valid cells
I_v_Holder = ones(1,nCells);
I_v_Holder(cIX) = 0;
I_v = find(I_v_Holder);

CR_dtr = CR_dtr_full(I_v,:);
STD = STD_full(I_v);
nCells = length(I_v);
CInfo = CInfo_full(I_v);

%% visualize valid cells, sorted by noise
[~,I] = sort(STD);
cIX = I;
gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
M = CR_dtr(I,1:1000);
BasicPlotMaps(cIX,gIX,M,CInfo_thr,stim,anat_yx,anat_yz,anat_zx);

%% Save new 'directload.mat'
temp = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load.mat']);
varList = {'CR_dtr','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn','perc_keep'}; 
% NOTE! not saving CR_raw in non-'_nodiscard.mat' again!
% should keep '_nodiscard.mat' so that this thresholding can be re-done
save(temp,varList{:},'-v7.3');



