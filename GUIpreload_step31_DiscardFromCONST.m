% (This is tailored for fishset>1)

clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_dir = GetCurrentDataDir();
save_dir = GetCurrentDataDir();
i_fish = 11;

%% load data
disp(['load fish ' num2str(i_fish)]);

filename = ['CONST_F' num2str(i_fish) '_fast_nodiscard.mat'];
[CellResp,const,dimCR] = LoadFileFromParts(data_dir,filename);

CellResp_full = CellResp;

%% Discard cells
isRankcells = true;

if isRankcells,
    % Method 1: discard __% noisy cells based on std of zscore of baseline
    % ! Manually set threshold:
    perc_keep = 50;
        
    % Compute std of dimest 10% of frames for each cell.
    tic
    prc = prctile(CellResp_full,10,2);
    nCells = size(CellResp_full,1);
    STD_full = zeros(nCells,1);
    for i = 1:nCells,
        ix = find(CellResp_full(i,:)<prc(i));
        STD_full(i) = std(CellResp_full(i,ix));
    end
    toc
    beep
    
    %% threshold and discard
    thr = prctile(STD_full,perc_keep);
    
    temp = find(STD_full>thr);
    [~,I] = sort(STD_full(temp));
    cIX = temp(I);
    
    %% visualize cells to discard
    % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
    gIX = round(cIX/1000)+1; % option: view anatomical z index
    M = CellResp_full(cIX,1:1000);
    stim = const.stim_full;
    CInfo = const.CInfo;
    anat_yx = const.anat_yx;
    anat_yz = const.anat_yz;
    anat_zx = const.anat_zx;
    BasicPlotMaps(cIX,gIX,M,CInfo,stim,anat_yx,anat_yz,anat_zx);
    
    % I_v: index of valid cells
    I_v_Holder = ones(1,nCells);
    I_v_Holder(cIX) = 0;
    I_v = find(I_v_Holder);
    
    CellResp = CellResp_full(I_v,:);
    STD = STD_full(I_v);
    nCells = length(I_v);
    CInfo_thr = CInfo(I_v);
    
    %% visualize valid cells, sorted by noise
    [~,I] = sort(STD);
    cIX = I;
    gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
    M = CellResp(I,1:1000);
    BasicPlotMaps(cIX,gIX,M,CInfo_thr,stim,anat_yx,anat_yz,anat_zx);
    
    %% update
    const.CInfo = CInfo(I_v);
else
    % Method 2: discard half undiscriminatedly (unranked)
    nCells = size(CellResp_full,1);
    cIX = (1:2:nCells)';
    
    %% visualize cells
    
    % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
    gIX = round(cIX/1000)+1; % option: view anatomical z index
    M = CellResp_full(cIX,1:1000);
    stim = const.stim_full;
    CInfo = const.CInfo;
    anat_yx = const.anat_yx;
    anat_yz = const.anat_yz;
    anat_zx = const.anat_zx;
    BasicPlotMaps(cIX,gIX,M,CInfo,stim,anat_yx,anat_yz,anat_zx);
    
    % update
    CellResp = CellResp_full(cIX,:);
    
    const.CInfo = CInfo(cIX);
    
end



%% Save new .mat - remember to rename properly!
% with partitioning of main data
newfishdir = fullfile(save_dir,['CONST_F' num2str(i_fish) '_fast_50noise.mat']);
dimCR = size(CellResp);
save(newfishdir,'const','dimCR','-v6');
SaveFileInPartsAppendv6(newfishdir,CellResp);


