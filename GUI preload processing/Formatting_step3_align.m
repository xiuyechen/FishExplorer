% Formatting_step3_align

clear all;close all;clc

%% set manually
M_stimset = GetFishStimset();
M_dir = GetFishDirectories();
save_masterdir = GetNestedDataDir();

%% 
range_fish = GetFishRange();

for i_fish = range_fish,
    disp(['i_fish = ', num2str(i_fish)]);

    %% load data
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(save_dir,'CoreInfo.mat'));
    load(fullfile(save_dir,'OptionalInfo.mat'));
    load(fullfile(save_dir,'AdditionalInfo.mat'));

    %% index processing
    if M_stimset(i_fish)==1, % fish 1-7
        % 'periods'
        M_period = {480,160,160,320,480,140,150,8,9,10,11,12,13,14,15,200}; %,{120,150,360}};% 8-15 are place holders
        periods = M_period{i_fish};
        
        % 'timelists'
        period = periods;        
        nrep = floor(size(frame_keys,1)/period)-1;
        shift = period;        
        timelists_raw = 1+shift:period*nrep+shift;
        
        % 'timelists_name'
        if i_fish<6,
            timelists_names = {'16 permut: B/W/phototaxis*2'};
        else % i_fish = 6 or 7
            timelists_names = {'phototaxis'};
        end
        
        % stimuluskey_raw
        stimuluskey_raw = round(frame_keys(:,17)');
        
    else % fish recorded with multiple protocols                      
        % 'stimset': struct that contains all info parsed from recorded params
        [stimset,stimuluskey_raw] = StimulusKeyPreprocessing(frame_keys,i_fish); % StimulusKeyPreprocessing(frame_keys,i_fish,'isplotting') to show plot
                
        % 'periods','timelists_names'
        periods = [stimset.period];
        timelists_names = {stimset.name};

        % very first 100 frames of experiment should be discarded (prone to artifacts)
        if stimset(1).starts(1)<100,
            if stimset(1).nReps(1)>1,
                stimset(1).starts(1) = stimset(1).starts(1)+stimset(1).period;
                stimset(1).nReps(1) = stimset(1).nReps(1)-1;
            elseif stimset(1).nReps(1)==1, % (n/a for current fish, improbable for future fish)
                stimset(1).ij(1,:) = [];
                stimset(1).rawstarts(1) = [];
                stimset(1).rawstops(1) = [];
                stimset(1).starts(1) = [];
                stimset(1).stops(1) = [];
                stimset(1).nReps(1) = [];
            end
        end
        
        % 'timelists' = list of time-frames numbers for a certain stimulus
        nTypes = length(stimset);
        timelists_raw = cell(1,nTypes); % time-lists, i.e. frame indices      
        for i_type = 1:nTypes,
            nSets = length(stimset(i_type).starts);
            IX = [];
            for i_set = 1:nSets,
                IX = horzcat(IX,stimset(i_type).starts(i_set):stimset(i_type).stops(i_set));
            end
            timelists_raw{i_type} = IX;
        end           
    end    
        
    %% Behavioral data: fictive motor recordings
    rows = [7,8,9,13,14];
    Behavior_raw = frame_keys(:,rows)';    
    Behavior_raw(2,:) = -Behavior_raw(2,:);
          
    %% save to nested directory
    save_dir_nested = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    
    filename = fullfile(save_dir_nested,'CoreInfo.mat');
    varList = {'timelists_raw','timelists_names','periods','stimuluskey_raw'}; % {'CellXYZ','anat_stack','fpsec'};
    save(filename,varList{:},'-append');
    
    filename = fullfile(save_dir_nested,'OptionalInfo.mat');
    varList = {'Behavior_raw'}; % {'numcell_full','CellXYZ_norm'};
    save(filename,varList{:},'-append');    
    
    if exist('stimset','var'),
        filename = fullfile(save_dir_nested,'AdditionalInfo.mat');
        save(filename,'stimset','-append'); % {'frame_keys','IX_inval_anat'};
    end    
end

%     % renameing
%     Behavior_raw = F;
%     stimuluskey_raw = stim_full_raw;
%     timelists = tlists_raw;
%     timelists_names = stimrangenames;

