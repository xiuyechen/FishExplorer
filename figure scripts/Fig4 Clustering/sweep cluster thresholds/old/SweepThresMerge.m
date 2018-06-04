data_masterdir = GetCurrentDataDir();

range_fish = 1:15;
% M_ClusGroup = [2,2,2,2];
% M_Cluster = [1,1,1,1];

% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2];

M_Thres = zeros(length(range_fish),3);
fract_dist = 0.8;

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
    LoadFullFish(hfig,i_fish);
    
    absIX = getappdata(hfig,'absIX');
    
    %% Cluster indexing
    i_ClusGroup = 1;%M_ClusGroup(i);
    i_Cluster = 1;%M_Cluster(i);
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    numK = Cluster(i_Cluster).numK;
    gIX = Cluster(i_Cluster).gIX;
    
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    setappdata(hfig,'cIX',cIX);
    toc
    
    %% load M_0;
    periods = getappdata(hfig,'periods');
    if length(periods)>1,
        setappdata(hfig,'stimrange',1:1);%length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    UpdateTimeIndex(hfig); % set M_0
    M_0 = getappdata(hfig,'M_0');
    
    M = M_0(cIX,:);
    
    %%
    M_Thres(i_fish,1) = FindMergeThresFromPDist(fract_dist,M);
    M_Thres(i_fish,2) = FindMergeThresFromPDist(fract_dist,M(:,1:round(end/2)));
    M_Thres(i_fish,3) = size(M,2);
    
    
    
figure;plot(M_Thres(:,3),M_Thres(:,1),'o')
    %%
%     f.UpdateIndices(hfig,cIX,gIX);
%     f.RefreshFigure(hfig);
end

