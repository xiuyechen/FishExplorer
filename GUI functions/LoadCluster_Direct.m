function [cIX,gIX,numK] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX)
if ~exist('absIX','var'),    
    data_masterdir = GetCurrentDataDir();
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
    absIX = h5read(hdf5_dir,'/absIX');
end
    
global VAR
ClusGroup = VAR(i_fish).ClusGroup{i_ClusGroup};
if length(ClusGroup)>=i_Cluster,
    numK = ClusGroup(i_Cluster).numK;
    gIX = ClusGroup(i_Cluster).gIX;
    cIX_abs = ClusGroup(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
else
    cIX = [];
    gIX = [];
    numK = [];
end
end