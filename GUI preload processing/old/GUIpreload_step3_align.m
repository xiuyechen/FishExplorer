%%%%%%%%%%%%%%%
%%
if true,
    clear all;close all;
    varList = {'CR_dtr','numcell_full','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
else
    % or without 'CR_dtr', keeping it from previous step:
    clearvars -except 'CR_dtr'; % 'CR_z_dtr';
    varList = {'numcell_full','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};
end

%%
M_dir = GetFishDirectories();
M_stimset = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2];

save_dir = GetCurrentDataDir();

dshift = 2; % differential shift between fictive and fluo, manually chosen

%% MANUAL
for i_fish = 1,
    disp(['i_fish = ', num2str(i_fish)]);
    tic
    %% load data
    load(fullfile(M_dir{i_fish},['Fish' num2str(i_fish) '_direct_load_full_v6.mat']),varList{:}); % '_direct_load_nodiscard.mat'
%     load(fullfile(datadir,'frame_turn.mat'),'frame_turn');
    
    if size(frame_turn,1)<size(frame_turn,2),
        frame_turn = frame_turn';
    end    
    
    %% Manual occasional frame corrections: ONLY RUN ONCE!!!!!!!!!
    
    if false,
        % add 1 frame at start: ------------- for which fish again??
        CR_dtr = horzcat(CR_dtr(:,1),CR_dtr);
        frame_turn = vertcat(frame_turn(1,:),frame_turn);
    end

    % correction of an error in Fish #1
    if i_fish == 1,
        CR_dtr = horzcat(CR_dtr(:,1:20),CR_dtr);
%         CR_z_dtr = horzcat(CR_z_dtr(:,1:20),CR_z_dtr);
    end
    
    %% index processing
    if M_stimset(i_fish)==1, % fish 1-7
        M_period = {480,160,160,320,480,140,150}; %,{120,150,360}};
        periods = M_period{i_fish};
        period = periods;
        
        nrep = floor(size(CR_dtr,2)/period)-1;

        shift = period; % (fixed 9/23/15; before = period+dshift)
        
        IX_all = 1+shift:period*nrep+shift;
        CellResp = CR_dtr(:,IX_all);
%         CellRespZ = CR_z_dtr(:,IX_all); % z-scored
        
        nCells = size(CR_dtr,1);
        CellRespAvr = mean(reshape(CR_dtr(:,IX_all+dshift),nCells,period,[]),3);
%         CellRespAvrZ = mean(reshape(CR_z_dtr(:,IX_all+dshift),nCells,period,[]),3);
        
        if i_fish<6,
            stimrangenames = {'16 permut: B/W/phototaxis*2'};
        else % i_fish = 6 or 7
            stimrangenames = {'phototaxis'};
        end
        %         stimrangenames = {'rep average','all reps','rep #1','rep #2','last rep'};
        ix_avr = 1:period;
        ix_all = 1:period*nrep;
        tlists = {ix_avr, ix_all};
        for i = [1,2,nrep],
            ix = 1+period*(i-1):period*i;
            tlists = [tlists, ix];
        end
        
        stim_full = frame_turn(IX_all,17)';
        stim_full = round(stim_full);
        
        stimAvr = stim_full(1:period);
        
    else % multiple protocols                
        % process raw stimulus code
        [stimset,stim_full_raw] = StimulusKeyPreprocessing(frame_turn,i_fish); % StimulusKeyPreprocessing(frame_turn,i_fish,'isplotting') to show plot
                
        %% variables to save later in struct
        periods = [stimset.period];
        stimrangenames = {stimset.name};

        %% find time-lists = list of time-frames numbers for a certain stimulus
        nTypes = length(stimset);
        tlists_raw = cell(1,nTypes+1); % time-lists, i.e. frame indices
        IX_all = [];
        for i_type = 1:nTypes,
            nSets = length(stimset(i_type).starts);
            IX = [];
            for i_set = 1:nSets,
                IX = horzcat(IX,stimset(i_type).starts(i_set):stimset(i_type).stops(i_set));
            end
            tlists_raw{i_type} = IX;
            IX_all = horzcat(IX_all,IX);
        end
        tlists_raw{nTypes+1} = IX_all;
%         tlists_raw{nTypes+2} = sort(IX_all);
        
        % This is the main data to store and load to GUI
        CellResp = CR_dtr(:,IX_all);
%         CellRespZ = CR_z_dtr(:,IX_all);
        stim_full = stim_full_raw(IX_all);
        
        %% get tlists corresponding to IX_all (tlists_raw ~ 1:totalframe#)
        tlists = cell(1,nTypes+1); % time-lists, i.e. frame indices
        for i_type = 1:nTypes+1,
            [~,IX] = ismember(tlists_raw{i_type},IX_all);
            tlists{i_type} = IX(find(IX));
        end
        
        %% find average
        numcell = size(CR_dtr,1);
        shift = -dshift; % circshift fluo left by 2 frames
        
        CellRespAvr = [];
        CellRespAvrZ = [];
        stimAvr = [];
        for i = 1:nTypes,
            M = circshift(CellResp(:,tlists{i}),shift,2);
            CellRespAvr = horzcat(CellRespAvr,mean(reshape(M,numcell,periods(i),[]),3));
%             CellRespAvrZ = horzcat(CellRespAvrZ,mean(reshape(M,numcell,periods(i),[]),3));
            % stim from stim_full
            m = stim_full(tlists{i});
            stimAvr = horzcat(stimAvr,m(1:periods(i)));
        end

    end
    
    %% prepare fictive data
    rows = [7,8,9,13,14];
    F = frame_turn(:,rows)';
    
    F(2,:) = -F(2,:);
    Fictive = F(:,IX_all);
    
    % find averages
    if M_stimset(i_fish)==1, % fish 1-7
        FictiveAvr = mean(reshape(Fictive,length(rows),period,[]),3);
    else
        FictiveAvr = [];
        for i = 1:nTypes,
            avr = mean(reshape(F(:,tlists_raw{i}),length(rows),periods(i),[]),3);
            FictiveAvr = horzcat(FictiveAvr,avr);
        end
    end
    
    % normalizations
    for i = 1:3,
        m = FictiveAvr(i,:);
        FictiveAvr(i,:) = (m-min(m))/(max(m)-min(m));
        m = Fictive(i,:);
        Fictive(i,:) = (m-min(m))/(max(m)-min(m));
    end
    m = FictiveAvr(4:5,:);
    FictiveAvr(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    m = Fictive(4:5,:);
    Fictive(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    
    
    %% compile CONST
    CellRespAvrZ = zscore(CellRespAvr')';
    
    data = [];
    names = {'numcell_full','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','periods','shift','dshift',...
        'CellRespAvr','CellRespAvrZ','Fictive','FictiveAvr','stim_full','stimAvr',... %'CellResp',
        'tlists','stimrangenames'};
    if M_stimset(i_fish) > 1,
        names = [names,{'tlists_raw','stimset'}];
    end
    
    for i = 1:length(names), % use loop to save variables into fields of CONST
        eval(['data.',names{i},'=', names{i},';']);
    end
    
    %% and save
    
    %%% new method with partitioning of main data
    newfishdir = fullfile(save_dir,['Data_F' num2str(i_fish) '_full.mat']);
    dimCR = size(CellResp);
    save(newfishdir,'data','dimCR','-v6');
    % custom function:
    SaveFileInPartsAppendv6(newfishdir,CellResp,'CellResp');
%     SaveFileInPartsAppendv6(newfishdir,CellRespZ,'CellRespZ');

    toc; beep
    
end


%% compute z-score
% disp('compute z-score')
% cell_resp_z_full = zscore(cell_resp_full')';

