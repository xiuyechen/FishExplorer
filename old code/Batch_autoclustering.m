% batch run full-clustering on all fish

% f.pushbutton_autoclus_Callback

global hm1;
hObject = hm1;

data_masterdir = GetCurrentDataDir();

range_fish = [5,6,7];
M_ClusGroup = [2,2,2,2];
M_Cluster = [1,1,1,1];

% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]; 

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);

    f.LoadFullFish(hfig,i_fish);
    absIX = getappdata(hfig,'absIX');
    
    %% Cluster indexing
    i_ClusGroup = M_ClusGroup(i);
    i_Cluster = M_Cluster(i);
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    numK = Cluster(i_Cluster).numK;
    gIX = Cluster(i_Cluster).gIX;
    
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    setappdata(hfig,'cIX',cIX);
    toc
            
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