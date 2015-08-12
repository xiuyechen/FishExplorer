%%%%%%%%%%%%%%%
%%
% clear all;close all;
% varList = {'CR_dtr','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn','periods'};

% % or
varList = {'CR_dtr','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn'};

%%
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

M_stimset = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2];

save_dir = 'C:\Janelia2014';

%% MANUAL
for i_fish = 10, %:8,    
    disp(num2str(i_fish));
    tic
    %% load data
    datadir = M_dir{i_fish};
    load(fullfile(datadir,['Fish' num2str(i_fish) '_direct_load.mat']),varList{:});
    
    %% ONLY RUN ONCE!!!!!!!!!
    % add 1 frame at start:
    CR_dtr = horzcat(CR_dtr(:,1),CR_dtr);

    frame_turn = vertcat(frame_turn(1,:),frame_turn);
    
    % correction of an error in Fish #1
    if i_fish == 1,
        CR_dtr = horzcat(CR_dtr(:,1:20),CR_dtr);
    end
    
    %% index processing
    if M_stimset(i_fish)==1, % fish 1-7, simple (single) protocol
        M_period = {480,160,160,320,480,280,300}; %,{120,150,360}};
        periods = M_period{i_fish};
        period = periods;
        
        nrep = floor(size(CR_dtr,2)/period)-1;
        dshift = 2;
        shift = period+dshift;
        
        IX_all = 1+shift:period*nrep+shift;
        CellResp = CR_dtr(:,IX_all); % old name: CRZt, 'Cell Responses Zscore trimmed'
        
        IX_avr = 1+shift:period+shift;
        nCells = size(CR_dtr,1);
        CellRespAvr = mean(reshape(CR_dtr(:,IX_all),nCells,period,[]),3); % 'Cell Response Average Zscore'
        
        datanames = {'rep average','all reps','rep #1','rep #2','last rep'};
        ix_avr = 1:period;
        ix_all = 1:period*nrep;
        tlists = {ix_avr, ix_all};
        for i = [1,2,nrep],
            ix = 1+period*(i-1):period*i;
            tlists = [tlists, ix];
        end
        
        stim_full = frame_turn(17,IX_all);
        stim_full = round(stim_full);
        
    else % multiple protocols
        % parse raw stimulus code
        [stimset,stim_full_raw] = StimulusKeyPreprocessing(frame_turn,i_fish);
        
        % find time-lists = list of time-frames numbers for a certain stimulus
        nTypes = length(stimset);
        tlists_raw = cell(1,nTypes+1); % time-lists, i.e. frame indices % each individual field, plus 'all'
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
        
        % This is the main data to store and load to GUI
        CellResp = CR_dtr(:,IX_all); % old name: CRZt, 'Cell Responses Zscore trimmed'
        stim_full = stim_full_raw(IX_all);
        
        %% get tlists corresponding to IX_all (tlists_raw ~ 1:totalframe#)
        tlists = cell(1,nTypes+1); % time-lists, i.e. frame indices
        for i_type = 1:nTypes+1,
            [~,tlists{i_type}] = intersect(IX_all,tlists_raw{i_type});
        end        
        
        %% variables to save later in struct % alternatively could compress/cleanup
        periods = [stimset.period];
        datanames = {'rep average','all reps',stimset.name};  
        
        dshift = 2;
        shift = -dshift; % circshift fluo left by 2 frames
        
        %% find average
%         CellRespAvr = []; % 'Cell Response Average Zscore'
%         stimtypelist = 1:nTypes;
%         numcell = size(CR_dtr,1);
%         for i = stimtypelist,
%             M = circshift(CR_dtr(:,tlists_raw{i}),shift,2);
%             CellRespAvr = horzcat(CellRespAvr,mean(reshape(M,numcell,periods(i),[]),3));
%         end                                
        
    end
    
    %% prepare fictive data
    rows = [7,8,9,13,14];
    F = frame_turn(:,rows)';
    
    F(2,:) = -F(2,:);
    Fc = F(:,IX_all);
        
    % find averages
    if M_stimset(i_fish)==1, % fish 1-7
        FcAvr = mean(reshape(Fc,length(rows),period,[]),3);
    else
        FcAvr = [];
        for i = 1:length(tlists)-1,
            avr = mean(reshape(Fc(:,tlists{i}),length(rows),periods(i),[]),3);
            FcAvr = horzcat(FcAvr,avr);
        end
    end
    
    % normalizations
    for i = 1:3,
        m = FcAvr(i,:);
        FcAvr(i,:) = (m-min(m))/(max(m)-min(m));
        m = Fc(i,:);
        Fc(i,:) = (m-min(m))/(max(m)-min(m));
    end
    m = FcAvr(4:5,:);
    FcAvr(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    m = Fc(4:5,:);
    Fc(4:5,:) =  (m-min(min(m)))/(max(max(m))-min(min(m)));
    
    
    %% compile CONST
    CONST = [];
    names = {'ave_stack','anat_yx','anat_yz','anat_zx','CInfo','periods','shift','CellResp',... %'CellRespAvr',
        'dshift','Fc','FcAvr','stim_full','stimset','tlists','datanames'};
    for i = 1:length(names), % use loop to save variables into fields of CONST
        eval(['CONST.',names{i},'=', names{i},';']);
    end
    
    %% and save
    %%% old method:
    %     temp = whos('CONST');
    %     if [temp.bytes]<2*10^9,
    %         save(['CONST_F' num2str(i_fish) '.mat'],'CONST');
    %         beep;
    %     else
    %         save(['CONST_F' num2str(i_fish) '.mat'],'CONST','-v7.3');
    %         beep;
    %     end
    
    %%% new method with partitioning of main data
    newfishdir = fullfile(save_dir,['CONST_F' num2str(i_fish) '_fast.mat']);
    const = CONST;
    const = rmfield(const,'CellResp');
    dimCR = size(CONST.CellResp);
    save(newfishdir,'const','dimCR','-v6');
    % custom function:
    SaveFileInPartsAppendv6(newfishdir,CellResp);

    toc; beep
    
end

    %% initialize VAR (once)
    

%% Initialize VARS % outdated? Clusgroup?

nCells = length(CONST.CInfo);

    cIX = (1:nCells)';
    i = 3;
    VAR(i_fish).Class(i).name = 'all processed';
    VAR(i_fish).Class(i).cIX = cIX;
    VAR(i_fish).Class(i).gIX = ones(length(cIX),1);
    VAR(i_fish).Class(i).numel = nCells;
    VAR(i_fish).Class(i).numK = 1;
    VAR(i_fish).Class(i).datatype = 'std';
    
    cIX = (1:100:nCells)';
    VAR(i_fish).ClusGroup{1,1}.name = 'test';
    VAR(i_fish).ClusGroup{1,1}.cIX = cIX;
    VAR(i_fish).ClusGroup{1,1}.gIX = ones(length(cIX),1);
    VAR(i_fish).ClusGroup{1,1}.numel = length(cIX);
    VAR(i_fish).ClusGroup{1,1}.numK = 1;

%     %%
%         cIX = (1:10:nCells)';
%     VAR(i_fish).ClusGroup{1,1}(12).name = '1/10 of all';
%     VAR(i_fish).ClusGroup{1,1}(12).cIX = cIX;
%     VAR(i_fish).ClusGroup{1,1}(12).gIX = ones(length(cIX),1);
%     VAR(i_fish).ClusGroup{1,1}(12).numel = length(cIX);
%     VAR(i_fish).ClusGroup{1,1}(12).numK = 1;

%     %% varience/std for reps for each cell
% %     if i_fish==2 || i_fish==3 || i_fish==6,
% %         period_real = period/2;
% %     else
% %         period_real = period;
% %     end
% %     nrep_real = floor((size(CR,2)-shift)/period_real);
% %     while period_real*nrep_real+shift>size(CR,2),
% %         nrep_real = nrep_real-1;
% %     end
% %     CRZ_3D = reshape(CRZ(:,1+shift:period_real*nrep_real+shift),nCells,period_real,[]);
% %     %% updated method, weighing both std between each rep and (summed with) std of 1st half & 2nd half of experiment - 1/8/15
% %     % CRZ = CONST.M_array.CellResp;
% %     % if i_fish==2 || i_fish==3 || i_fish==6,
% %     %     period_real = CONST.M_array.period/2;
% %     % else
% %     %     period_real = CONST.M_array.period;
% %     % end
% %     % CRZ_3D = reshape(CRZ,size(CRZ,1),period_real,[]);
% %     % divide = round(size(CRZ_3D,3)/2);
% %     % CRZ_std1 = std(CRZ_3D(:,:,1:divide),0,3);
% %     % CRZ_std2 = std(CRZ_3D(:,:,divide+1:end),0,3);
% %     % temp1 = mean(CRZ_std1,2);
% %     % temp2 = mean(CRZ_std2,2);
% %     %
% %     % temp12 = horzcat(temp1,temp2);
% %     % temp = mean(temp12,2)+std(temp12,0,2);
% %     % [~,I] = sort(temp);
% %     % M = temp(I);
% %     % figure;plot(M)
% %     %
% %     % figure;imagesc(CRZ(I,:))
% %     %
% %     % nCells = size(CRZ,1);
% %     
% %     %% find low variance / stimulus-locked cells
% %     CRZ_std = std(CRZ_3D,0,3);
% %     temp = mean(CRZ_std,2);
% %     
% %     % find mean-std thres: 0.5
% %     [~,I] = sort(temp);
% %     M = temp(I);
% %     figure;plot(M)
% %     %%
% %     i_last = length(VAR(i_fish).Class);
% %     M_perc = [0.025,0.1,0.3];
% %     for j = 1:length(M_perc);
% %         thres = M(round(nCells*M_perc(j)));
% %         cIX = find(temp<thres);
% %         i = j+i_last;
% %         VAR(i_fish).Class(i).round = 0;
% %         VAR(i_fish).Class(i).name = ['perc < ' num2str(M_perc(j)*100) '%'];
% %         %     VAR(i_fish).Class(i).notes = ['perc < ' num2str(M_perc(j)*100) '%'];
% %         VAR(i_fish).Class(i).cIX = cIX;
% %         VAR(i_fish).Class(i).gIX = ones(length(cIX),1);
% %         VAR(i_fish).Class(i).numel = length(cIX);
% %         VAR(i_fish).Class(i).numK = 1;
% %         VAR(i_fish).Class(i).datatype = 'std';
% %     end
% %     
% %     %% shift CR?
% %     % shift = 161;
% %     % nrep=floor(size(CR,2)/period)-1;
% %     %
% %     % skiplist=[];
% %     % IX_rep=setdiff(1:nrep, skiplist);
% %     % IX=zeros(period*length(IX_rep),1);
% %     %
% %     % for i=1:length(IX_rep)
% %     %     IX(period*(i-1)+1:period*i)=period*(IX_rep(i)-1)+1+shift:period*(IX_rep(i))+shift;
% %     % end
% %     %
% %     % % cell_resp=cell_resp(:,1:nrep*period);
% %     % CRA=mean(reshape(CR(:,IX),[nCells,period,length(IX_rep)]),3);
% %     % CRAZ = zscore(CRA')';
% %     
% 
