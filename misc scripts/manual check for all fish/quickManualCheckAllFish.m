% Check all fish!

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

range_fish = 1:18;
data_masterdir = GetCurrentDataDir();
M_stimrange = GetStimRange();%'O');
M_fishset = GetFishStimset();

TF = zeros(length(range_fish),1);
for i_fish = range_fish,    
    LoadFullFish(hfig,i_fish,0);
%     absIX = getappdata(hfig,'absIX');
    numcell_full = getappdata(hfig,'numcell_full');
    TF(i_fish,1) = numcell_full;%length(absIX) ==numcell_full;   
end