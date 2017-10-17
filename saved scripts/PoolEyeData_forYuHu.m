
clear all;clc;
%  close all; 
outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = [1:12,14:18];% stricter range: [1:8,11,12,14:17]; % fish 12 left is not great
ClusterIDs = [12,1]; % ABN (seed)

S = [];

for i = 1:length(range_fish)
    i_fish = range_fish(i);
    [cIX_eye,gIX_eye,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    S(i).cIX_eye = cIX_eye;
    S(i).cIX_eye = cIX_eye;
    % should be similar to Behavior_full; BehaviorAvr
    S(i).Eye_full_trace = getappdata(hfig,'Eye_full_motorseed');
    S(i).Eye_avr_trace = getappdata(hfig,'EyeAvr_motorseed');
    S(i).absIX = getappdata(hfig,'absIX');
    S(i).M = M;
end

save(fullfile(outputDir,'Eye_data.mat'),'S','range_fish');