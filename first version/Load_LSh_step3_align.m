%%%%%%%%%%%%%%%
%%
clear all;close all;
%%
M_dir = {'F:\Janelia2014\Fish1_16states_30frames';
    'F:\Janelia2014\Fish2_20140714_2_4_16states_10frames';
    'F:\Janelia2014\Fish3_20140715_1_1_16_states_10frames';
    'F:\Janelia2014\Fish4_20140721_1_8_16states_20frames';
    'F:\Janelia2014\Fish5_20140722_1_2_16states_30frames';
    'F:\Janelia2014\Fish6_20140722_1_1_3states_30,40frames';
    'F:\Janelia2014\Fish7_20140722_2_3_3states_30,50frames';
    'F:\Janelia2014\Fish8_20141222_2_2_7d_PT_3OMR_shock_lowcut';
    'F:\Janelia2014\Fish9_20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356'};

% global VAR;
% load('VAR_current.mat');

% load('CONSTS.mat');
% VAR = [];
% CONST = [];

M_set = [1, 1, 1, 1, 1, 1, 1, 2];

%% MANUAL
for i_fish = 9, %:8,    
    disp(num2str(i_fish));
    %% load data
    datadir = M_dir{i_fish};
    load(fullfile(datadir,['Fish' num2str(i_fish) '_direct_load.mat']));
    %     varList = {'CR_raw','CR_dtr','nCells','CInfo','anat_yx','anat_yz','ave_stack','fpsec','frame_turn','periods'};
    
    %% ONLY RUN ONCE!!!!!!!!!
    % add 1 frame at start:
    CR_dtr = horzcat(CR_dtr(:,1),CR_dtr);
    frame_turn = horzcat(frame_turn(:,1),frame_turn);

    % correction of an error in Fish #1
    if i_fish == 1,
        CR_dtr = horzcat(CR_dtr(:,1:20),CR_dtr);
    end
    
    %% index processing
    if M_set(i_fish)==1, % fish 1-7
        period = periods{1};
        nrep = floor(size(CR_dtr,2)/period)-1;
        dshift = 2;
        shift = period+dshift;
        
        IX_all = 1+shift:period*nrep+shift;
        CRZt = CR_dtr(:,IX_all); % 'Cell Responses Zscore trimmed'
        
        IX_avr = 1+shift:period+shift;
        nCells = size(CR_dtr,1);
        CRAZ = mean(reshape(CR_dtr(:,IX_all),nCells,period,[]),3);
        
        datanames = {'rep average','all reps','rep #1','rep #2','last rep'};
        ix_avr = 1:period;
        ix_all = 1:period*nrep;
        tlists = {ix_avr, ix_all};
        for i = [1,2,nrep],
            ix = 1+period*(i-1):period*i;
            tlists = [tlists, ix];
        end
    else % fish 8...
        % (manually segmented/extracted from frame_turn(17,:))               
        dshift = 2;
        shift = -dshift; % circshift fluo left by 2 frames
        
        period_pt = 120; % = periods{1} % from direct load...
        period_omr = 150; % = periods{2}
        period_spt = 360; % = periods{3}
        
        %% unshifted
        IX_pt1 = 1:1080; % first session
        IX_pt2 = 5581:6660; % second session
        IX_pt = horzcat(IX_pt1,IX_pt2);
        
        IX_pt2_circ = circshift(IX_pt2,60,2); % shift of 60 determined with corr
        IX_pt_circ = horzcat(IX_pt1,IX_pt2_circ);
        
        IX_omr = horzcat(1081:2430,6661:8010);
        IX_spt = horzcat(2431:2790,5221:5580);
        IX_ptomr = union(IX_pt,IX_omr);
        
        numrep = length(IX_pt)/period_pt; % = length(IX_omr)/period_omr = 18
        IX_ptomr_circ = []; % interlaced for easy averaging with period = period_pt+period_omr!
        for i = 1:numrep,
            IX_ptomr_circ = [IX_ptomr_circ,IX_pt_circ((i-1)*period_pt+1:i*period_pt)];
            IX_ptomr_circ = [IX_ptomr_circ,IX_omr((i-1)*period_omr+1:i*period_omr)];
        end
        IX_all = union(IX_ptomr, IX_spt);

        %% shift, to get CRZt and CRAZ
        IX_pt1_ = circshift(1:1080,shift,2); % first session 
        IX_pt2_ = circshift(5581:6660,shift,2); % second session
        IX_pt_ = horzcat(IX_pt1_,IX_pt2_);   
        
        IX_pt2_circ_ = circshift(IX_pt2_,60,2); % shift of 60 determined with corr                   
        IX_pt_circ_ = horzcat(IX_pt1_,IX_pt2_circ_);
        
        IX_omr1_ = circshift(1081:2430,shift,2); % first session 
        IX_omr2_ = circshift(6661:8010,shift,2); % second session 
        IX_omr_ = horzcat(IX_omr1_,IX_omr2_);
        IX_spt1_ = circshift(2431:2790,shift,2); % first session
        IX_spt2_ = circshift(5221:5580,shift,2); % second session
        IX_spt_ = horzcat(IX_spt1_,IX_spt2_);
        IX_ptomr_ = union(IX_pt_,IX_omr_);
        
        numrep = length(IX_pt_)/period_pt; % = length(IX_omr)/period_omr = 18
        IX_ptomr_circ_ = []; % interlaced for easy averaging with period = period_pt+period_omr!
        for i = 1:numrep,
            IX_ptomr_circ_ = [IX_ptomr_circ_,IX_pt_circ_((i-1)*period_pt+1:i*period_pt)];
            IX_ptomr_circ_ = [IX_ptomr_circ_,IX_omr_((i-1)*period_omr+1:i*period_omr)];
        end        
        IX_all_ = union(IX_ptomr_, IX_spt_);             
        CRZt = CR_dtr(:,IX_all_); % 'Cell Responses Zscore trimmed'
        
        nCells = size(CR_dtr,1);
        % find averages of different stimuli, and splice together
        a1 = mean(reshape(CR_dtr(:,IX_pt_circ_),nCells,period_pt,[]),3);
        a2 = mean(reshape(CR_dtr(:,IX_omr_),nCells,period_omr,[]),3);
%         a3 = mean(reshape(CR_dtr(:,IX_spt),nCells,period_spt,[]),3);
        CRAZ = horzcat(a1,a2);%,a3); % average of only PT and OMR for now       

        %% time lists (from unshifted)
        datanames = {'rep average','all reps','phototaxis','OMR grating','spontaneous (black)','phototaxis+OMR'};
        [~,ix_pt_circ] = ismember(IX_pt_circ,IX_all);
        [~,ix_omr] = ismember(IX_omr,IX_all);
        [~,ix_spt] = ismember(IX_spt,IX_all);
        [~,ix_ptomr_circ] = ismember(IX_ptomr_circ,IX_all);
        %         ix_avr = horzcat(ix_pt_(1:period_pt),ix_omr(1:period_omr),ix_spt(1:period_spt));
        ix_avr_ptomr = horzcat(ix_pt_circ(1:period_pt),ix_omr(1:period_omr));
        tlists = {ix_avr_ptomr,1:size(CRZt,2),ix_pt_circ,ix_omr,ix_spt,ix_ptomr_circ};
        
    end
    
    %% prepare fictive data
    rows = [7,8,9,13,14];
    F = frame_turn(rows,:);
    
    F(2,:) = -F(2,:);
    Fc = F(:,IX_all);
    
    % find averages
    if M_set(i_fish)==1, % fish 1-7
        FcAvr = mean(reshape(Fc,length(rows),period,[]),3);
    else
        a1 = mean(reshape(F(:,IX_pt_circ),length(rows),period_pt,[]),3);
        a2 = mean(reshape(F(:,IX_omr),length(rows),period_omr,[]),3);
%         a3 = mean(reshape(F(:,IX_spt),length(rows),period_spt,[]),3);
        FcAvr = horzcat(a1,a2);%,a3);
    end
    
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
        
    photostate = frame_turn(17,IX_all);     
    if M_set(i_fish)==1,
        photostate = round(photostate);
    else
        M1 = [2,3,4,5,12,13,14,15,99]; % original values in photostate
        M2 = [3,1,3,2,4 ,10,11,12, 0]; % standardized values
        ps = photostate; % copy
        % for transitional values, "round" to target photostates
        ps_rounded = interp1(M1,M1,photostate,'nearest');
        % swap
        for i = 1:length(M1),
            photostate(ps_rounded==M1(i)) = M2(i);
        end
    end
    
    %% compile CONST
    CONST = [];
    names = {'ave_stack','anat_yx','anat_yz','anat_zx','CInfo','periods','shift','CRAZ','CRZt','dshift','FcAvr','Fc','photostate','tlists','datanames'};
    for i = 1:length(names), % use loop to save variables into fields of CONST
        eval(['CONST.',names{i},'=', names{i},';']);
    end
    
    %% initialize VAR (once)
    
    %%
    save(['CONST_F' num2str(i_fish) '.mat'],'CONST');
    beep;
    %% OR
    save(['CONST_F' num2str(i_fish) '.mat'],'CONST','-v7.3');
    beep;
end



%% Initialize VARS % outdated? Clusgroup?

nCells = length(CONST.CInfo);

    cIX = (1:nCells)';
    i = 1;
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

    %%
        cIX = (1:10:nCells)';
    VAR(i_fish).ClusGroup{1,1}(12).name = '1/10 of all';
    VAR(i_fish).ClusGroup{1,1}(12).cIX = cIX;
    VAR(i_fish).ClusGroup{1,1}(12).gIX = ones(length(cIX),1);
    VAR(i_fish).ClusGroup{1,1}(12).numel = length(cIX);
    VAR(i_fish).ClusGroup{1,1}(12).numK = 1;

    %% varience/std for reps for each cell
%     if i_fish==2 || i_fish==3 || i_fish==6,
%         period_real = period/2;
%     else
%         period_real = period;
%     end
%     nrep_real = floor((size(CR,2)-shift)/period_real);
%     while period_real*nrep_real+shift>size(CR,2),
%         nrep_real = nrep_real-1;
%     end
%     CRZ_3D = reshape(CRZ(:,1+shift:period_real*nrep_real+shift),nCells,period_real,[]);
%     %% updated method, weighing both std between each rep and (summed with) std of 1st half & 2nd half of experiment - 1/8/15
%     % CRZ = CONST.M_array.CRZt;
%     % if i_fish==2 || i_fish==3 || i_fish==6,
%     %     period_real = CONST.M_array.period/2;
%     % else
%     %     period_real = CONST.M_array.period;
%     % end
%     % CRZ_3D = reshape(CRZ,size(CRZ,1),period_real,[]);
%     % divide = round(size(CRZ_3D,3)/2);
%     % CRZ_std1 = std(CRZ_3D(:,:,1:divide),0,3);
%     % CRZ_std2 = std(CRZ_3D(:,:,divide+1:end),0,3);
%     % temp1 = mean(CRZ_std1,2);
%     % temp2 = mean(CRZ_std2,2);
%     %
%     % temp12 = horzcat(temp1,temp2);
%     % temp = mean(temp12,2)+std(temp12,0,2);
%     % [~,I] = sort(temp);
%     % M = temp(I);
%     % figure;plot(M)
%     %
%     % figure;imagesc(CRZ(I,:))
%     %
%     % nCells = size(CRZ,1);
%     
%     %% find low variance / stimulus-locked cells
%     CRZ_std = std(CRZ_3D,0,3);
%     temp = mean(CRZ_std,2);
%     
%     % find mean-std thres: 0.5
%     [~,I] = sort(temp);
%     M = temp(I);
%     figure;plot(M)
%     %%
%     i_last = length(VAR(i_fish).Class);
%     M_perc = [0.025,0.1,0.3];
%     for j = 1:length(M_perc);
%         thres = M(round(nCells*M_perc(j)));
%         cIX = find(temp<thres);
%         i = j+i_last;
%         VAR(i_fish).Class(i).round = 0;
%         VAR(i_fish).Class(i).name = ['perc < ' num2str(M_perc(j)*100) '%'];
%         %     VAR(i_fish).Class(i).notes = ['perc < ' num2str(M_perc(j)*100) '%'];
%         VAR(i_fish).Class(i).cIX = cIX;
%         VAR(i_fish).Class(i).gIX = ones(length(cIX),1);
%         VAR(i_fish).Class(i).numel = length(cIX);
%         VAR(i_fish).Class(i).numK = 1;
%         VAR(i_fish).Class(i).datatype = 'std';
%     end
%     
%     %% shift CR?
%     % shift = 161;
%     % nrep=floor(size(CR,2)/period)-1;
%     %
%     % skiplist=[];
%     % IX_rep=setdiff(1:nrep, skiplist);
%     % IX=zeros(period*length(IX_rep),1);
%     %
%     % for i=1:length(IX_rep)
%     %     IX(period*(i-1)+1:period*i)=period*(IX_rep(i)-1)+1+shift:period*(IX_rep(i))+shift;
%     % end
%     %
%     % % cell_resp=cell_resp(:,1:nrep*period);
%     % CRA=mean(reshape(CR(:,IX),[nCells,period,length(IX_rep)]),3);
%     % CRAZ = zscore(CRA')';
%     

