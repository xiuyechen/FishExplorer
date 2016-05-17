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
    
    %% load M_0;
    periods = getappdata(hfig,'periods');
    if length(periods)>1,
        setappdata(hfig,'stimrange',1:3);%length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    UpdateTimeIndex(hfig); % set M_0
    M_0 = getappdata(hfig,'M_0');
    
    %% Cluster indexing
    i_ClusGroup = 1;%M_ClusGroup(i);
    i_Cluster = 3;%M_Cluster(i);
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    numK = Cluster(i_Cluster).numK;
    gIX = Cluster(i_Cluster).gIX;
    
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
    setappdata(hfig,'cIX',cIX);
    toc
    
    M = M_0(cIX,:);
    
    cIX_20 = cIX;
    gIX_20 = gIX;
    
    %%
    range_numK = 10:10:60;
    range_thres = 0.4:0.025:0.7;
    Stats = zeros(length(range_numK),length(range_thres));
    for i_numK = 1:length(range_numK),
        cIX = cIX_20;
        gIX = gIX_20;
        % divide
        numK2 = range_numK(i_numK);
        gIX = KmeansSubdivide(numK2,gIX,M_0);
        
        % foxel regression
        Reg = FindCentroid_Direct(gIX,M);
        [cIX,gIX,~] = AllCentroidRegression_direct(M_0,0.75,Reg);
        gIX = SqueezeGroupIX(gIX);
        
        % size threshold
        U = unique(gIX);
        numU = length(U);
        for i=1:numU,
            if length(find(gIX==U(i)))<thres_minsize/2,
                cIX(gIX==U(i)) = [];
                gIX(gIX==U(i)) = [];
            end
        end
        [gIX,numU] = SqueezeGroupIX(gIX);
        disp(numU);
        
        % histogram?
        
        C = FindCentroid_Direct(gIX,M_0(cIX,:));
        
        Ccorr = 1-pdist(C,'correlation');
        Ccorr_sq = squareform(Ccorr);
        hist(Ccorr,-1:0.05:1)
        
        for i_thres = 1:length(range_thres),
            
            IX = find(Ccorr>range_thres(i_thres));
            [IX_ctr,~] = ind2sub(size(Ccorr_sq),IX);
            list = unique(IX_ctr);
            
            Count = 0;
            for i_ctr = 1:length(IX_ctr),
                temp = find(gIX==IX_ctr(i_ctr));
                Count = Count + length(temp);
            end
            
            Stats(i_numK,i_thres) = length(IX_ctr)/2/length(Ccorr_sq); % avr connection per cluster
        end
    end
    %%
    figure; hold on
    M_title = {'# cell','# clus','% of clus'};
    for j = 1:length(range_numK),
        subplot(length(range_numK),1,j)
        plot(range_thres,Stats(j,:))
        title(['numK=',num2str(20*range_numK(j))]);
    end
    %%
    thres_minsize = 10;
    thres_merge = 0.8;
    thres_reg = thres_merge;%0.7;
    thres_cap = 0.7;
    
    %% 1.2. Regression with the centroid of each cluster
    disp('regression with all clusters');
    tic
    Reg = FindCentroid_Direct(gIX,M);
    [cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg,Reg);
    gIX = SqueezeGroupIX(gIX);
    toc
    % clusgroupID = 1;
    % SaveCluster_Direct(cIX,gIX,absIX,i_fish,'k20x20_reg',clusgroupID);
    % size thresholding
    U = unique(gIX);
    numU = length(U);
    for i=1:numU,
        if length(find(gIX==U(i)))<thres_minsize/2,
            cIX(gIX==U(i)) = [];
            gIX(gIX==U(i)) = [];
        end
    end
    [gIX,numU] = SqueezeGroupIX(gIX);
    disp(numU);
    %% Find Seed
    tic
    [cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,thres_reg,cIX,gIX,M_0);
    toc
    
    %     isWkmeans = 1;
    %     AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans)
    
    %%
    f.UpdateIndices(hfig,cIX,gIX);
    f.RefreshFigure(hfig);
end



%         U = unique(gIX);
%         H = zeros(length(U),4);
%         for i = 1:length(U),
%             IX = find(gIX==U(i));
%             Dist = pdist(M_0(cIX(IX),:),'corr');
%             cdist = corr(C(i,:)',M_0(cIX(IX),:)');
%             if ~isempty(Dist),
%                 H(i,1) = 1-mean(Dist);
%                 H(i,2) = 1-max(Dist);
%                 H(i,3) = min(cdist);
%             else
%                 H(i,:) = NaN;
%             end
%         end
%         %
%         figure; hold on
%         subplot(141);hist(H(:,1),-0.5:0.05:1)
%         subplot(142);hist(H(:,2),-0.5:0.05:1)
%         subplot(143);hist(H(:,3),-0.5:0.05:1)