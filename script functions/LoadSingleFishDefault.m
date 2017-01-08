function [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,h0)

setappdata(h0,'i_fish',i_fish);

%% load fish data
isFullData = true;
LoadFullFish(h0,i_fish,isFullData);

%% load the auto-clustering indices
i_ClusGroup = 6; % master-thres 0.7
i_Cluster = 1;

[cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);

setappdata(h0,'cIX',cIX);
setappdata(h0,'gIX',gIX);

%% load cluster time-series data
% set params
setappdata(h0,'isAvr',0);
setappdata(h0,'isRefAnat',1);
% setappdata(h0,'isZsore',0); % this is the default

% load
M_stimrange = GetStimRange();
stimrange = M_stimrange{i_fish};
timelists = getappdata(h0,'timelists');
fishset = getappdata(h0,'fishset');

tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
[M,behavior,stim] = GetTimeIndexedData_Default_Direct(h0,cIX,tIX);

setappdata(h0,'M',M);
setappdata(h0,'behavior',behavior);
setappdata(h0,'stim',stim);
setappdata(h0,'tIX',tIX);
setappdata(h0,'stimrange',stimrange);

M_0 = GetTimeIndexedData(h0,'isAllCells');
setappdata(h0,'M_0',M_0);
end