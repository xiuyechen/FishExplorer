% InitializationforGUI_step2_cellselection

clear all;close all;clc

save_masterdir = GetCurrentDataDir();

%%
range_fish = GetFishRange();

for i_fish = range_fish,
    %% load data
    tic
    disp(['load fish ' num2str(i_fish) '...']);    
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    % load 'data'
    load(fullfile(save_dir,'data_full.mat'),'data'); % struct with many fields
    data_full = data;
    names = fieldnames(data_full); % cell of strings
    for i = 1:length(names),
        eval([names{i} ' = data_full.' names{i} ';']);
    end
    % load time series (hdf5 file)
    CellRespZ_full = h5read(fullfile(save_dir,'TimeSeries.h5'),'/CellRespZ');
    toc
    
    %% generate absolute cell index 
    % (this indexing is stable; any later cell selection would be a subset of this)
    absIX_full = (1:numcell_full)';
    absIX = setdiff(absIX_full,IX_inval_anat);
    
    %% Discard cells
    % isRankcells = true;
    % %%
    % if isRankcells,
    % Method 1: discard __% noisy cells based on std of zscore of baseline
    
    % ! Manually set threshold:
    perc_keep = 50;
    
    % Compute std of dimest 10% of frames for each cell.
    disp('compute std');
    tic
    prc = prctile(CellRespZ_full,10,2);
    STD_full = zeros(length(absIX),1);
    for i = 1:length(absIX),
        ix = find(CellRespZ_full(i,:)<prc(i));
        STD_full(i) = std(CellRespZ_full(i,ix));
    end
    toc
    beep
    
    %% threshold and discard
    thr = prctile(STD_full,perc_keep);
    
%     if true,
        temp = find(STD_full>thr);
        [~,I] = sort(STD_full(temp));
        IX_inval_rank = absIX(temp(I));
        [~,cIX] = ismember(IX_inval_rank,absIX);        
        gIX = round(cIX/1000)+1; % option: view anatomical z index
        % gIX = (round((1:length(cIX))/1000)+1)'; % option: view sorted index
%     else
%         % contrast discard vs. non-discard
%         temp1 = find(STD_full>thr);
%         [~,I1] = sort(STD_full(temp1));
%         temp2 = find(STD_full<thr);
%         [~,I2] = sort(STD_full(temp2));
%         cIX = [temp1(I1);temp2(I2)];
%         gIX = [ones(length(I1),1);2*ones(length(I2),1)];
%     end
    %% visualize cells to discard
%     M = CellResp_full(cIX,1:1000);
%     BasicPlotMaps(cIX,gIX,M,CellXYZ,stim_full,anat_yx,anat_yz);
    
    %% visualize valid cells, sorted by noise
    %     [~,I] = sort(STD);
    %     cIX = I;
    %     gIX = (round((1:length(cIX))/1000)+1)'; % sorted index
    %     M = CellResp(I,1:1000);
    % %     BasicPlotMaps(cIX,gIX,M,CellXYZ,stim,anat_yx,anat_yz);
    
    %% update: trim cell indices
%     ClusGroup(clusID).cIX_abs = absIX(cIX);
%     [~,cIX] = ismember(cIX_abs,absIX);
%     IX_inval_rank = absIX(cIX);
    
    absIX_inval = union(IX_inval_anat,IX_inval_rank);
    absIX_half = setdiff(absIX_full,absIX_inval);
    [~,cIX_half] = ismember(absIX_half,absIX);
    
    % load other sets
    CellResp_full = h5read(fullfile(save_dir,'TimeSeries.h5'),'/CellResp');
    CellRespAvr_full = h5read(fullfile(save_dir,'TimeSeries.h5'),'/CellRespAvr');
    CellRespAvrZ_full = h5read(fullfile(save_dir,'TimeSeries.h5'),'/CellRespAvrZ');
    
    % trim all
    CellResp = CellResp_full(cIX_half,:);
    CellRespZ = CellRespZ_full(cIX_half,:);
    CellRespAvr = CellRespAvr_full(cIX_half,:);
    CellRespAvrZ = CellRespAvrZ_full(cIX_half,:);
    
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
    
    %%
    filename = fullfile(save_dir,'TimeSeries_half.h5');
    h5create(filename,'/CellResp',size(CellResp),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellResp',CellResp);
    
    h5create(filename,'/CellRespZ',size(CellRespZ),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespZ',CellRespZ);
    
    h5create(filename,'/CellRespAvr',size(CellRespAvr),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespAvr',CellRespAvr);
    
    h5create(filename,'/CellRespAvrZ',size(CellRespAvrZ),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespAvrZ',CellRespAvrZ);
    
    h5create(filename,'/absIX',size(absIX_half),'Datatype','single');
    h5write(filename,'/absIX',absIX_half);
    
end

