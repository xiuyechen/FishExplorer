% Generate VAR to save all clusters

clear all;close all;clc

data_masterdir = GetCurrentDataDir();
data_masterdir;

range_fish = GetFishRange();

%% Initialize VAR
VAR = [];
% load(fullfile(save_dir,'VAR_new.mat'),'VAR');
%%
% for i_fish = range_fish,
%     %% load data
%     disp(['load fish ' num2str(i_fish) '...']);
%     data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
%     load(fullfile(data_dir,'data_full.mat'),'data'); % struct with many fields
% 
%     hdf5_dir = fullfile(data_dir,'TimeSeries_half.h5');
%     absIX_half = h5read(hdf5_dir,'/absIX');
%         
%     hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
%     absIX = h5read(hdf5_dir,'/absIX');
%     
%     %%
%     i_ClusGroup = 1;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'selection';
%     
%     i_Cluster = 1;
%     cIX_abs = absIX_half(1:10:end);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '1/10 of 50%_rank';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
%     
%     %
%     i_ClusGroup = 2;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'init';   
%     
%     i_Cluster = 1;
%     cIX_abs = absIX_half;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '50%_rank cells';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
%     
%     i_Cluster = 2;
%     cIX_abs = absIX; 
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'all valid cells';
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
%     VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;
%     
%     %
%     i_ClusGroup = 3;
%     VAR(i_fish).ClusGroupName{i_ClusGroup} = 'AutoClusters';  
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
%     
% end
% 
% %%
% save(fullfile(save_dir,'VAR_new.mat'),'VAR');
% disp('saved updated ''VAR''');
% 
