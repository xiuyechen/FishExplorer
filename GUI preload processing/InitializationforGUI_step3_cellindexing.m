% InitializationforGUI_step3_cellindexing

% initialize VAR
clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_masterdir = GetCurrentDataDir();
save_dir = data_masterdir;

range_fish = GetFishRange();

%% Initialize VAR
% VAR = [];
load(fullfile(save_dir,'VAR_new.mat'),'VAR');
%%
for i_fish = range_fish,
    %% load data
    disp(['load fish ' num2str(i_fish) '...']);
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(data_dir,'data_full.mat'),'data'); % struct with many fields

    hdf5_dir = fullfile(data_dir,'TimeSeries_half.h5');
    absIX_half = h5read(hdf5_dir,'/absIX');
        
    hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
    absIX = h5read(hdf5_dir,'/absIX');
    
    %%
    i_ClusGroup = 1;
    VAR(i_fish).ClusGroupName{i_ClusGroup} = 'Selection';
    
    i_Cluster = 1;
    cIX_abs = absIX_half(1:10:end);
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '1/10 of 50%_rank';
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
    
    %
    i_ClusGroup = 2;
    VAR(i_fish).ClusGroupName{i_ClusGroup} = 'Init';   
    
    i_Cluster = 1;
    cIX_abs = absIX;
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'all valid cells';
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
    
    i_Cluster = 2;
    cIX_abs = absIX_half;
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '50%_rank cells';
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
    VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
    
    %
%     i_ClusGroup = 3;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'AutoClusters_raw';  
%     
%     i_Cluster = 1;
%     cIX_abs = absIX_half;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'temp';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
%     
%     %
%     i_ClusGroup = 4;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'figures';  
%     
%     i_Cluster = 1;
%     cIX_abs = absIX_half(1:10:end);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'temp';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
%     
%     %   
%     i_ClusGroup = 5;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'other';  
%     
%     i_Cluster = 1;
%     cIX_abs = absIX_half(1:10:end);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'temp';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
    
end

%%
save(fullfile(save_dir,'VAR_new.mat'),'VAR');
disp('saved updated ''VAR''');

%%
% varience/std for reps for each cell
%% updated method, weighing both std between each rep and (summed with) std of 1st half & 2nd half of experiment - 1/8/15
% CRZ = CONST.M_array.CellResp;
% if i_fish==2 || i_fish==3 || i_fish==6,
%     period_real = CONST.M_array.period/2;
% else
%     period_real = CONST.M_array.period;
% end
% CRZ_3D = reshape(CRZ,size(CRZ,1),period_real,[]);
% divide = round(size(CRZ_3D,3)/2);
% CRZ_std1 = std(CRZ_3D(:,:,1:divide),0,3);
% CRZ_std2 = std(CRZ_3D(:,:,divide+1:end),0,3);
% temp1 = mean(CRZ_std1,2);
% temp2 = mean(CRZ_std2,2);
%
% temp12 = horzcat(temp1,temp2);
% temp = mean(temp12,2)+std(temp12,0,2);
% [~,I] = sort(temp);
% M = temp(I);
% figure;plot(M)
%
% figure;imagesc(CRZ(I,:))
%
% nCells = size(CRZ,1);

%% find low variance / stimulus-locked cells
% CRZ_std = std(CRZ_3D,0,3);
% temp = mean(CRZ_std,2);
% 
% % find mean-std thres: 0.5
% [~,I] = sort(temp);
% M = temp(I);
% figure;plot(M)
% %%
% i_last = length(VAR(i_fish).Class);
% M_perc = [0.025,0.1,0.3];
% for j = 1:length(M_perc);
%     thres = M(round(nCells*M_perc(j)));
%     cIX = find(temp<thres);
%     i = j+i_last;
%     VAR(i_fish).Class(i).round = 0;
%     VAR(i_fish).Class(i).name = ['perc < ' num2str(M_perc(j)*100) '%'];
%     %     VAR(i_fish).Class(i).notes = ['perc < ' num2str(M_perc(j)*100) '%'];
%     VAR(i_fish).Class(i).cIX = cIX;
%     VAR(i_fish).Class(i).gIX = ones(length(cIX),1);
%     VAR(i_fish).Class(i).numel = length(cIX);
%     VAR(i_fish).Class(i).numK = 1;
%     VAR(i_fish).Class(i).datatype = 'std';
% end


