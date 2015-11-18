
clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_dir = GetCurrentDataDir();
save_dir = GetCurrentDataDir();
i_fish = 1;

%% load data
disp(['load fish ' num2str(i_fish)]);

filename = ['Data_F' num2str(i_fish) '_full.mat'];
load(fullfile(data_dir,filename),'data');
[CellResp_full,data_full,dimCR] = LoadFileFromParts(data_dir,filename,'CellResp');
% [CellRespZ_full] = LoadFileFromParts(data_dir,filename,'CellRespZ');

% load anat outliers, or move code here...
 load('C:\Janelia2014\Fish1_16states_30frames\Fish1_extrainfo_anat.mat','IX_inval_anat');
 
 numcell_full = length(data_full.CellXYZ);
 absIX_full = (1:numcell_full)';
 absIX = setdiff(absIX_full,IX_inval_anat); 
 
%% Discard cells
% isRankcells = true;
% %%
% if isRankcells,    
    % Method 1: discard __% noisy cells based on std of zscore of baseline
            
    % ! Manually set threshold:
    perc_keep = 50;
    
    disp('compute z-score');
    tic
    CellRespZ_full = zscore(CellResp_full')';
    toc
    % Compute std of dimest 10% of frames for each cell.
    disp('compute std');
    tic
    prc = prctile(CellRespZ_full(absIX,:),10,2);
    STD_full = zeros(length(absIX),1);
    for i = 1:length(absIX),
        ix = find(CellRespZ_full(absIX(i),:)<prc(i));
        STD_full(i) = std(CellRespZ_full(absIX(i),ix));
    end
    toc
    beep
    
    %% threshold and discard
    thr = prctile(STD_full,perc_keep);
    
    if true,
        temp = find(STD_full>thr);
        [~,I] = sort(STD_full(temp));
        IX_inval_rank = absIX(temp(I));
        cIX = IX_inval_rank;
        gIX = round(cIX/1000)+1; % option: view anatomical z index
        % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
    else
        % contrast discard vs. non-discard
        temp1 = find(STD_full>thr);
        [~,I1] = sort(STD_full(temp1));
        temp2 = find(STD_full<thr);
        [~,I2] = sort(STD_full(temp2));
        cIX = [temp1(I1);temp2(I2)];
        gIX = [ones(length(I1),1);2*ones(length(I2),1)];
    end
    %% visualize cells to discard
    M = CellResp_full(cIX,1:1000);
    stim = data_full.stim_full;
    CellXYZ = data_full.CellXYZ;
    anat_yx = data_full.anat_yx;
    anat_yz = data_full.anat_yz;
    BasicPlotMaps(cIX,gIX,M,CellXYZ,stim,anat_yx,anat_yz);

    %% visualize valid cells, sorted by noise
%     [~,I] = sort(STD);
%     cIX = I;
%     gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
%     M = CellResp(I,1:1000);
% %     BasicPlotMaps(cIX,gIX,M,CellXYZ,stim,anat_yx,anat_yz);
    
    %% update
    absIX_inval = union(IX_inval_anat,IX_inval_rank);
    absIX = setdiff(absIX_full,absIX_inval);
    
    CellResp = CellResp_full(absIX,:);
    CellRespZ = CellRespZ_full(absIX,:);
    
    data = data_full;
    data.CellXYZ = CellXYZ(absIX,:);    
    data.CellRespAvr = data_full.CellRespAvr(absIX,:);
    data.CellRespAvrZ = data_full.CellRespAvrZ(absIX,:);
    data.absIX = absIX;    
    data.IX_inval_anat = IX_inval_anat;
    data.IX_inval_rank = IX_inval_rank;
    
% else
%     % Method 2: discard half undiscriminatedly (unranked)
%     nCells = size(CellResp_full,1);
%     cIX = (1:2:nCells)';
%     
%     %% visualize cells
%     
%     % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
%     gIX = round(cIX/1000)+1; % option: view anatomical z index
%     M = CellResp_full(cIX,1:1000);
%     stim = data_full.stim_full;
%     CellXYZ = data_full.CellXYZ;
%     anat_yx = data_full.anat_yx;
%     anat_yz = data_full.anat_yz;
%     anat_zx = data_full.anat_zx;
%     BasicPlotMaps(cIX,gIX,M,CInfo,stim,anat_yx,anat_yz,anat_zx);
% end

%% Save new .mat - remember to name properly!
% with partitioning of main data
newfishdir = fullfile(save_dir,['Data_F' num2str(i_fish) '.mat']);
dimCR = size(CellResp);
save(newfishdir,'data','dimCR','-v6');
SaveFileInPartsAppendv6(newfishdir,CellResp,'CellResp');
SaveFileInPartsAppendv6(newfishdir,CellRespZ,'CellRespZ');


