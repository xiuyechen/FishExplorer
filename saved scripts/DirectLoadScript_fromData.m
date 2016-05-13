data_masterdir = GetCurrentDataDir();

range_fish = [5,6,7];
M_ClusGroup = [2,2,2,2];
M_Cluster = [1,1,1,1];

% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]; 

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);

    LoadFullFish(hfig,i_fish);
    
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
    UpdateTimeIndex(hfig); % set M_0
    M_0 = getappdata(hfig,'M_0');
    
    %%
    f.UpdateIndices(hfig,cIX,gIX,numK);
    
    
    %% 
    isWkmeans = 1;
    AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans)
end