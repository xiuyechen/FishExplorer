% InitializationforGUI_step1
% Processing nested data: one-time initialization for GUI

clear all; close all

behav_fluo_shift = -2; % shift correction between fictive and fluo, manually chosen


%% Load data
data_masterdir = GetNestedDataDir();
range_fish = GetFishRange();

for i_fish = range_fish,
    
    disp(['fish ' num2str(i_fish) ' loading...']);
    tic
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    
    TimeSeries = h5read(fullfile(data_dir,'TimeSeries.h5'),'/TimeSeries');
    
    filename = fullfile(data_dir,'CoreInfo.mat');
    load(filename);% {'periods','timelists','timelists_names','stimuluskey_raw','CellXYZ','anat_stack','fpsec'};
    
    filename = fullfile(data_dir,'OptionalInfo.mat');
    load(filename); % {'Behavior_raw','numcell_full','CellXYZ_norm','IX_inval'};
    
    filename = fullfile(data_dir,'AdditionalInfo.mat');
    load(filename); % {'frame_keys','IX_inval_anat'} and if exist, 'stimset';
    toc
    
    %% generate absolute cell index
    % (this indexing is stable; any later cell selection would be a subset of this)
    absIX_full = (1:numcell_full)';
    absIX = setdiff(absIX_full,IX_inval_anat);
    
    %% Anatomy stack projections
    % x-y view
    im = max(anat_stack,[],3);
    out=imNormalize99(im);
    anat_yx = repmat(out,[1 1 3]);
    
    % y-z view
    im = squeeze(max(anat_stack,[],2));
    out=imNormalize99(im);
    anat_yz = repmat(out,[1 1 3]);
    
    % x-z view
    im = squeeze(max(anat_stack,[],1));
    out=imNormalize99(im);
    out = flipud(out'); %%%% empirically necessary...
    anat_zx = repmat(out,[1 1 3]);
    
%     dimv_yx = size(anat_yx);
%     dimv_yz = size(anat_yz);
%     dimv_zx = size(anat_zx);
    
    %% get full-length time-index from timelists
    if length(timelists_names)==1, % fish with single stimulus set
        IX_all_raw = timelists_raw;
        % get relative time-index (relative to full-length list), ~ frame number
        timelists = {1:length(IX_all_raw)};
        
        % stim & behavior
        stim_full = stimuluskey_raw(IX_all_raw);
        Behavior_full = Behavior_raw(:,IX_all_raw);
        
        % trimming & shift
        CellResp_full = TimeSeries(:,IX_all_raw);
        CellResp_full = circshift(CellResp_full,behav_fluo_shift,2);
        
        % averages
        period = periods;
        numcell = size(CellResp_full,1);
        CellRespAvr_full = mean(reshape(CellResp_full,numcell,period,[]),3);
        stimAvr = stim_full(1:period);
        BehaviorAvr = mean(reshape(Behavior_full,size(Behavior_full,1),period,[]),3);
        
    else % fish with multiple stimulus sets
        
        temp = cat(2, timelists_raw{:});
        IX_all_raw = sort(temp);
        
        % get relative time-index (relative to full-length list), ~ frame number
        nTypes = length(timelists_names);
        timelists = cell(1,nTypes);
        for i_type = 1:nTypes,
            [~,IX] = ismember(timelists_raw{i_type},IX_all_raw);
            timelists{i_type} = IX(find(IX));
        end
        
        % stim & behavior
        stim_full = stimuluskey_raw(IX_all_raw);
        Behavior_full = Behavior_raw(:,IX_all_raw);
        
        % trimming & shift
        CellResp_full = TimeSeries(:,IX_all_raw);
        CellResp_full = circshift(CellResp_full,behav_fluo_shift,2); % using circshift as a shortcut for padding on one end
        
        % averages
        numcell = size(CellResp_full,1);
        CellRespAvr_full = [];
        stimAvr = [];
        BehaviorAvr = [];
        for i = 1:nTypes,
            M = CellResp_full(:,timelists{i});
            if mod(numel(M),numcell*periods(i))==0,
                avr = mean(reshape(M,numcell,periods(i),[]),3);
            else % for patterns with dummy periods where the stimulus is constant, like 'spont'
                nrep = floor(numel(M)/(numcell*periods(i)));
                M2 = M(:,1:nrep*periods(i));
                avr = mean(reshape(M2,numcell,periods(i),[]),3);
            end
            CellRespAvr_full = horzcat(CellRespAvr_full,avr); %#ok<AGROW>
            
            m = stim_full(timelists{i});
            stimAvr = horzcat(stimAvr,m(1:periods(i))); %#ok<AGROW>
            
            M_behav = Behavior_full(:,timelists{i});
            if mod(numel(M_behav),size(M_behav,1)*periods(i))==0,            
                avr = mean(reshape(M_behav,size(M_behav,1),periods(i),[]),3);
            else % for patterns with dummy periods where the stimulus is constant, like 'spont'
                nrep = floor(numel(M_behav)/(size(M_behav,1)*periods(i)));
                M2_behav = M_behav(:,1:nrep*periods(i));
                avr = mean(reshape(M2_behav,size(M_behav,1),periods(i),[]),3);
            end
            BehaviorAvr = horzcat(BehaviorAvr,avr); %#ok<AGROW>
        end
    end
    
    %% Behavior normalizations
    for i = 1:3,
        m = BehaviorAvr(i,:);
        BehaviorAvr(i,:) = (m-min(m))/(max(m)-min(m));
        m = Behavior_full(i,:);
        Behavior_full(i,:) = (m-min(m))/(max(m)-min(m));
    end
    m = BehaviorAvr(4:5,:);
    BehaviorAvr(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    m = Behavior_full(4:5,:);
    Behavior_full(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    
    
    %% compute z-score
    disp('compute z-score...')
    tic
    CellRespZ_full = zscore(CellResp_full')';
    CellRespAvrZ_full = zscore(CellRespAvr_full')';
    toc
    
    %% save
    disp('saving to files...')
    tic
    % directory to save data formatted for distribution:
    save_masterdir = GetCurrentDataDir();
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % save time-series
    filename = fullfile(save_dir,'TimeSeries.h5');
    h5create(filename,'/CellResp',size(CellResp_full(absIX,:)),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellResp',CellResp_full(absIX,:));
    
    h5create(filename,'/CellRespZ',size(CellRespZ_full(absIX,:)),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespZ',CellRespZ_full(absIX,:));
    
    h5create(filename,'/CellRespAvr',size(CellRespAvr_full(absIX,:)),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespAvr',CellRespAvr_full(absIX,:));
    
    h5create(filename,'/CellRespAvrZ',size(CellRespAvrZ_full(absIX,:)),'Datatype','single','ChunkSize',[1000 100]);
    h5write(filename,'/CellRespAvrZ',CellRespAvrZ_full(absIX,:));
    
    h5create(filename,'/absIX',size(absIX),'Datatype','single');
    h5write(filename,'/absIX',absIX);
    toc
    %%
    data = [];
    names = {'periods','timelists_names','stimuluskey_raw','CellXYZ','anat_stack','fpsec',...
        'Behavior_raw','numcell_full','CellXYZ_norm','IX_inval_anat',...
        'anat_yx','anat_yz','anat_zx',...
        'timelists','stim_full','stimAvr','Behavior_full','BehaviorAvr'};
%         'absIX'}; % absIX now stored in hdf5
    if length(timelists_names)>1, % M_stimset(i_fish) > 1,
        names = [names,{'stimset'}]; %#ok<AGROW>
    end
    
    for i = 1:length(names), % use loop to save variables into fields of 'data'
        eval(['data.',names{i},'=', names{i},';']);
    end
    
    save(fullfile(save_dir,'data_full.mat'),'data');
    
    toc
    
end
