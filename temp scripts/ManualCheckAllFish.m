% Check all fish!
range_fish = 1:18;
data_masterdir = GetCurrentDataDir();
M_stimrange = GetStimRange();%'O');
M_fishset = GetFishStimset();

TF = zeros(length(range_fish),2);
for i_fish = range_fish,    
    LoadFullFish(hfig,i_fish,-1);
    absIX = getappdata(hfig,'absIX');
    numcell_full = getappdata(hfig,'numcell_full');
    TF(i_fish,1) = numcell_full;%length(absIX) ==numcell_full;
    
    stimrange = M_stimrange{i_fish};
    fishset = M_fishset(i_fish);
    timelists = getappdata(hfig,'timelists');
    tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
    TF(i_fish,2) = length(tIX);
    
    [cIX,gIX] = LoadCluster_Direct(i_fish,7,1,absIX);
    TF(i_fish,3) = length(unique(gIX));
    
%     disp(['absIX:',num2str(length(absIX))]);
%     disp(['numcell_full:',num2str(numcell_full)]);
end