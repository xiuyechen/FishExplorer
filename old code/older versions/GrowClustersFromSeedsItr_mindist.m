function [cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,cIX,gIX,M_0)
%disp('find ROIs')

% Set params

% may be determined by hist of distances between all cells (in a sample)
% in correlation distance, i.e. 1 - corr.coeff
% thres_merge = 0.4;
% thres_cap = 0.5;

% thres_minsize = 10; % cell number in final clusters

%%
gIX = SqueezeGroupIX(gIX);
M = M_0(cIX,:);
C = FindCentroid_Direct(gIX,M); % import functional supervoxels ('foxel')
nFoxels = size(C,1);

%% Calculate correlation distance between all cluster-centroids
temp = pdist(C,'correlation');
Dist = squareform(temp);
for i = 1:nFoxels,
    Dist(i,i) = NaN;
end

%% ITERATION:
% initialize
ROI = [];
ROI(nFoxels).fxlist = []; % not capping
ROI(nFoxels).numcell = [];
ROIcount = 0;

% analogy ~ mask to store new ROI
Mask = zeros(nFoxels,1); 

%%

while true, % loop through all qualified seeds           
    % look for next seed (pair of foxels)
    [~,ix] = min(Dist(:));
    [I,J] = ind2sub(size(Dist),ix);
    % and add seed position to mask
    Mask(I) = 1;
    Mask(J) = 1;
    
    if Dist(I,J)>(1-thres_merge),
        break;
    end
    
    IX_list = [find(gIX==I);find(gIX==J)];
    M_core = M(IX_list,:);
    newCore = mean(M_core,1);
    
    while true % grow ROI: expand mask and examine next neighbors
        IX_rest = find(Mask==0);

        % find next closest foxel
        D = 1-corr(newCore',C(IX_rest,:)'); % distance matrix between new core and potential foxels
        if isempty(D),
            break;
        end
        [a,ix] = min(D);
        ix_cand = IX_rest(ix);
        IX_list_cand = [IX_list;find(gIX==ix_cand)];
        M_core_cand = M(IX_list_cand,:);
        newCore_cand = mean(M_core_cand,1);
        coredist = 1-corr(C(ix_cand,:)',newCore_cand');
        
        if a<(1-thres_merge) && coredist<(1-thres_cap),
            Mask(ix_cand) = 1;
            IX_list = IX_list_cand;
            newCore = newCore_cand;
        else
            break; % finished expanding this seed
        end
    end
    
    % save this ROI if bigger than thres
    numcell = length(IX_list);
    if numcell >= thres_minsize,
        ROIcount = ROIcount+1;
        ROI(ROIcount).fxlist = find(Mask);
        ROI(ROIcount).numcell = numcell;
    end
    
    % update/reset
    list = find(Mask);
    Dist(list,:) = NaN;
    Dist(:,list) = NaN;
    C(list,:) = NaN;
    Mask = zeros(nFoxels,1);
    
    if isempty(find(isnan(Dist)==0)),
        break;
    end
end

ROI(ROIcount+1:end) = [];

%% Merge accordingly
for i = 1:ROIcount,
    fxlist = ROI(i).fxlist;
    for j = 2:length(fxlist),
        gIX(gIX==fxlist(j)) = fxlist(1);
    end
end
[gIX,numU] = SqueezeGroupIX(gIX);

end

%% using average instead of min dist
% find next closest foxel: candidate foxel
