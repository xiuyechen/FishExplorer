% batch template
global hm1;
hObject = hm1;

data_masterdir = GetCurrentDataDir();

range_fish =  8:13;
% M_ClusGroup = [2,2,2,2];
% M_Cluster = [1,1,1,1];
const_ClusGroup = 2;
const_Cluster = 1;
% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]; 

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);

    f.LoadFullFish(hfig,i_fish);
    absIX = getappdata(hfig,'absIX');
    
    %% Cluster indexing
    i_ClusGroup = const_ClusGroup;% M_ClusGroup(i);
    i_Cluster = const_Cluster;% M_Cluster(i);
    ClusGroup = VAR(i_fish).ClusGroup{i_ClusGroup};
    numK = ClusGroup(i_Cluster).numK;
    gIX = ClusGroup(i_Cluster).gIX;    
    cIX_abs = ClusGroup(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    
    setappdata(hfig,'cIX',cIX);
            
    %%
    periods = getappdata(hfig,'periods');
    if length(periods)>1,
        setappdata(hfig,'stimrange',1:length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    UpdateTimeIndex(hfig);

    %%
    f.UpdateIndices(hfig,cIX,gIX,numK);
    
    %%
    tic
    pushbutton_autoclus_Callback(hObject,f,i_fish);
    toc
end