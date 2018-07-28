function [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,clusParams,isMakeFoxels,masterthres)
% AUTOCLUSTERING
% [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,[],isMakeFoxels,masterthres);
% Pipeline to automatically cluster all cells,
%
% INPUTS
% cIX: indices of currently selected cells; the algorithm starts with this selection. 
% gIX: grouping indices of currently selected cells (e.g. if k-means has been performed already, see 'isWkmeans')
% M_0: matrix of all neural activity (row ~ cells, col ~ time points)
% cIX_reg: scope for regression (correlation of traces),which indexes from the full array of cells in 'M_0')
% isWkmeans: flag for whether starting with performing k-means (1 or true) or using precomputed k-means clusters (0 or false)
% clusParams: (option) specify detailed parameters
%       clusParams = struct('merge',thres_merge,'cap',thres_cap,'reg1',thres_reg,...
%             'reg2',thres_reg2,'minSize',thres_minsize,'k1',numK1);
% isMakeFoxels: start with computing foxels (default true) instead of using precomputed foxels
% masterthres: (option) this replaces 4 thresholds in clusParams (merge, cap, reg1, reg2) with the same number (masterthres, default = 0.7)


%% Detailed description
%{
This algorithm was custom developed to suit this larval zebrafish dataset.
We outline the algorithm below:
1.	Divide all cells into “functional voxels” (~10 cells each)
    a.	Perform k-means clustering on all cells (k=20)
    b.	Perform k-means clustering on outputs of (a) (k = ~400)
    c.	Discard any cells whose correlation with the voxel average activity is less than $THRESH
    d.	Discard any voxels with fewer than 5 cells
2.	Merge voxels into clusters based on density in functional space
    a.	for each pair of voxels ij (starting from most correlated): 
    if the correlation between voxel i and j is greater than $THRESH, 
    and the correlation between the the voxel j and the centroid (average) 
    of the cluster containing voxel i is greater than $THRESH:
    then group voxel j in the same cluster as i. 
    b.	discard any clusters with fewer than 10 cells
3.	Clean up clusters using regression to cluster centroids
    a.	for each cell k:
    if the correlation to the closest cluster’s centroid is greater than $THRESH:
    include cell k in that cluster
    b.	Discard any clusters with fewer than 10 cells
4.	Iterate merge and cleanup steps 
    a.	Perform step 2 and 3 once more, using clusters as input voxels. 

This clustering algorithm can either be applied to all cells in the brain 
or a chosen subset of interest, and the correlation threshold determining 
clustering stringency ($THRESH) can be adjusted to trade-off completeness 
and accuracy (see Results and Fig. S4-1d). For most analysis in the text, 
the value of $THRESH was 0.7.

A whole-brain single-cell resolution data (~100,000 cells, ~5000 time frames) 
can be clustered on a standard desktop computer on the timescale of minutes.
The algorithm is sensitive enough to detect even very weak functional
clustering patterns in the data, including artefactual signals resulting 
from the scanning laser itself (Fig. S4-1i, these artefactual clusters were 
thereafter excluded from analysis).
%}

%% set params
% default:
if ~exist('isWkmeans','var')
    isWkmeans = false;
end
if ~exist('isMakeFoxels','var')
    isMakeFoxels = true;
end
if ~exist('masterthres','var')
    masterthres = 0.7;
end

% if isAutoclusWithAllCells,
%     cIX_reg = (1:size(M_0,1))';
% else
%     cIX_reg = cIX;
% end

thres_reg2 = masterthres;
thres_reg = masterthres;
thres_merge = masterthres;
thres_cap = masterthres;
thres_minsize = 10;
numK1 = 20;

isShowProgress = 1;

% optional override:
if exist('clusParams','var')
    if ~isempty(clusParams)       
        thres_merge = clusParams.merge; % correlation coeff
        thres_cap = clusParams.cap; % correlation coeff
        thres_reg2 = clusParams.reg2; % correlation coeff
        thres_minsize = clusParams.minSize; % number of cells
        
        isShowProgress = 0;
    else
        clusParams = struct('merge',thres_merge,'cap',thres_cap,'reg1',thres_reg,...
            'reg2',thres_reg2,'minSize',thres_minsize,'k1',numK1);
    end
else
    clusParams = struct('merge',thres_merge,'cap',thres_cap,'reg1',thres_reg,...
            'reg2',thres_reg2,'minSize',thres_minsize,'k1',numK1);
end

%% 1. Obtain foxels (functional voxels) 
% (density-based, effectively down-sampling)

autoClusStart = tic;

if isMakeFoxels % unless previously saved    
    [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,cIX_reg,isWkmeans,clusParams);%absIX,i_fish);
end

%% 2. Merge foxels to obtain fROI's (functional ROI's) 
% agglomerative clustering, performed in 'GrowClustersFromSeedsItr'
% (perform 2 rounds, inspired by iterative processes)

if isShowProgress
    disp('Merge foxels (round 1)');
end
grow1Start = tic;
M = M_0(cIX,:);
gIX = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,gIX,M);
grow1Time = toc(grow1Start);

% Regression with the centroid of each cluster (round 1)
reg1Start = tic;
Reg = FindCentroid_Direct(gIX,M_0(cIX,:));
[cIX,gIX,nMerge1] = AllCentroidRegression_SizeThres_direct(M_0(cIX_reg,:),cIX_reg,thres_reg2,Reg,thres_minsize/2);
reg1Time = toc(reg1Start);

if isShowProgress
    disp('Merge foxels (round 2)');
end

grow2Start = tic;
M = M_0(cIX,:);
gIX = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,gIX,M);
grow2Time = toc(grow2Start);

reg2Start = tic;
Reg = FindCentroid_Direct(gIX,M_0(cIX,:));
[cIX,gIX,nMerge2] = AllCentroidRegression_SizeThres_direct(M_0(cIX_reg,:),cIX_reg,thres_reg2,Reg,thres_minsize);
reg2Time = toc(reg2Start);

% clusgroupID = 2;
% SaveCluster_Direct(cIX,gIX,absIX,i_fish,'afterGrowth',clusgroupID);

%% clean-up and rank hierarchically for display
U = unique(gIX);
numU = length(U);
for i=1:numU
    if length(find(gIX==U(i)))<thres_minsize
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end

C = FindCentroid_Direct(gIX,M_0(cIX,:));
gIX = HierClus_Direct(C,gIX);
[gIX,numROI] = SqueezeGroupIX(gIX);

%% display timing
autoClusTime = toc(autoClusStart);

disp(['(Size) nMerge1:' num2str(nMerge1) ...
    ' nMerge2:' num2str(nMerge2) ...
    ' nROI:' num2str(numROI) ]);

disp(['(Time) grow1:' num2str(grow1Time)...
    ' reg1:' num2str(reg1Time) ...
    ' grow2:' num2str(grow2Time)...
    ' reg2:' num2str(reg2Time) ...
    ' total time:' num2str(autoClusTime) ' sec']);

end