function [cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,thres_reg,cIX,gIX,M_0)
disp('find ROIs')

% Set params

% may be determined by hist of distances between all cells (in a sample)
% in correlation distance, i.e. 1 - corr.coeff
% thres_merge = 0.4;
% thres_cap = 0.5;

% thres_minsize = 10; % cell number in final clusters

%%
gIX = SqueezeGroupIX(gIX);
M = M_0(cIX,:);
C = FindCentroid_Direct(gIX,M); % kmeans 20x20, 'supervoxels'
nPixels = size(C,1);

%% Calculate correlation distance between all cluster-centroids
temp = pdist(C,'correlation');
Dist = squareform(temp);
for i = 1:nPixels,
    Dist(i,i) = NaN;
end

%% ITERATION:
% initialize
ROI = [];
ROI(nPixels).pxlist = []; % not capping
ROI(nPixels).numcell = [];
ROIcount = 0;

% analogy ~ mask to store new ROI
BW = zeros(nPixels,1); 

%%
tic
while true, % loop through all qualified seeds           
    % look for next seed
    [~,ix] = min(Dist(:));
    [I,J] = ind2sub(size(Dist),ix);
    % and add seed position to mask
    BW(I) = 1;
    BW(J) = 1;
    
    if Dist(I,J)>(1-thres_merge),
        break;
    end
    
    while true % grow ROI: expand mask and examine next neighbors
        % find new core based on cells in clusters
        list = find(BW);
        IX = [];
        for i_gIX = 1:length(list),
            IX = [IX;find(gIX==list(i_gIX))];
        end
        M_core = M(IX,:);
        try
            [~,newCore] = kmeans(M_core,1,'distance','correlation');
        catch
            
        end
        
        % find next closest cluster
        D = 1-corr(newCore',C'); % distance matrix between new core and the original cluster means
        IX1 = find(BW==0);
        [a,IX2] = min(D(IX1));
        
        BW_pretend = BW;
        BW_pretend(IX1(IX2)) = 1;
        IX = find(BW_pretend);
        maxdist = max(pdist(C(IX,:),'correlation'));
        
        if a<(1-thres_merge) && maxdist<(1-thres_cap),
            BW(IX1(IX2)) = 1;
        else
            break; % finished expanding this seed
        end
    end
    
    % save this ROI if bigger than thres
    list = find(BW);
    IX = [];
    for i_gIX = 1:length(list),
        IX = [IX;find(gIX==list(i_gIX))]; %#ok<AGROW>
    end
    numcell = length(IX);
    if numcell >= thres_minsize,
        ROIcount = ROIcount+1;
        ROI(ROIcount).pxlist = find(BW);
        ROI(ROIcount).numcell = numcell;
    end
    
    % update/reset
    Dist(list,:) = NaN;
    Dist(:,list) = NaN;
    C(list,:) = NaN;
    BW = zeros(nPixels,1);
end

ROI(ROIcount+1:end) = [];
toc
%% Merge accordingly
for i = 1:ROIcount,
    pxlist = ROI(i).pxlist;
    for j = 2:length(pxlist),
        gIX(gIX==pxlist(j)) = pxlist(1);
    end
end
[gIX,numU] = SqueezeGroupIX(gIX);
disp(numU);

%% Regression with the centroid of each cluster, round 2
% disp('auto-reg');
% Reg = FindCentroid_Direct(gIX,M);
% [cIX,gIX] = AllCentroidRegression_direct(M_0,thres_reg,Reg);

%% size threshold
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX,numU] = SqueezeGroupIX(gIX);
disp(numU);
end

