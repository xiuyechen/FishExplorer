% save all centroids into single mat file
range_fish = 2:3;
%
% AllCentroids = cell(1,14);
i_ClusGroup = 3;
i_Cluster = 1;

for i_fish = range_fish;
    
    f.LoadFullFish(hfig,i_fish);
    
    periods = getappdata(hfig,'periods');
    if periods>1,
        setappdata(hfig,'stimrange',1:length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    
    setappdata(hfig,'cIX',1);
    setappdata(hfig,'gIX',1);
    UpdateTimeIndex(hfig);
    
    %% Cluster indexing
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    
    absIX = getappdata(hfig,'absIX');
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    %     cIX = cIX(find(cIX));
    gIX = Cluster(i_Cluster).gIX;
    numK = Cluster(i_Cluster).numK;
    %     setappdata(hfig,'cIX',cIX);
    %     setappdata(hfig,'gIX',gIX);
    %     setappdata(hfig,'numK',numK);
    f.UpdateIndices(hfig,cIX,gIX,numK);
    
    %% %     C = FindCentroid(hfig);
    
    M = getappdata(hfig,'M');    
    U = unique(gIX);
    numU = length(U);
    C = zeros(numU,size(M,2));
    D = zeros(numU,1);
    for i=1:numU,
        IX = find(gIX == U(i));
        if length(IX)==1,
            C(i,:) = M(IX,:);
            D(i) = 1;
        else
            M_s = M(IX,:);
            [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
            C(i,:) = C1;
            D(i) = mean(D1);
        end
    end
    
    %%
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    U = unique(gIX);
    numU = length(U);
    XYZn = cell(numU,1);
    for i = 1:numU,
        IX = find(gIX == U(i));
        XYZn{i} = CellXYZ_norm(cIX_abs(IX),:);
    end
    
    AllCentroids{i_fish}.Centroids = C;
    AllCentroids{i_fish}.XYZn = XYZn;
    AllCentroids{i_fish}.stim = getappdata(hfig,'stim');
    AllCentroids{i_fish}.behavior = getappdata(hfig,'behavior');
    
end

datadir = GetCurrentDataDir();
save(fullfile(datadir,'AllCentroids.mat'),'AllCentroids');