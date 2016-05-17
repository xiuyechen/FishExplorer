function [M_0,M,N,T,hfig] = getFullFishData(hfig,i_fish,VAR)

LoadFullFish(hfig,i_fish,0); %load 50% rank for now
    absIX = getappdata(hfig,'absIX');
    
     %% Cluster indexing
    i_ClusGroup = 2;%M_ClusGroup(i);
    i_Cluster = 1;%M_Cluster(i);
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    numK = Cluster(i_Cluster).numK;
    gIX = Cluster(i_Cluster).gIX;
    
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    setappdata(hfig,'cIX',cIX);
    toc
            
    %
    periods = getappdata(hfig,'periods');
    if length(periods)>1,
        setappdata(hfig,'stimrange',1:length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    UpdateTimeIndex(hfig); % set M_0
    
    
    M_0 = getappdata(hfig,'M_0');
    M = M_0(cIX,:);
    
    N = size(M_0,1);
    T = size(M_0,2);