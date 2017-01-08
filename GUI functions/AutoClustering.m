function [cIX,gIX] = AutoClustering(cIX,gIX,M_0,cIX_reg,isWkmeans,clusParams,isMakeFoxels,masterthres)%,isAutoclusWithAllCells)
% [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels,masterthres,isAutoclusWithAllCells)
% automatically cluster all cells, starting with the currently selected
% cells

%% set params
% default:
if ~exist('isWkmeans','var'),
    isWkmeans = false;
end
if ~exist('isMakeFoxels','var'),
    isMakeFoxels = true;
end
if ~exist('masterthres','var'),
    masterthres = 0.7;
end
% if ~exist('isAutoclusWithAllCells','var'),
%     isAutoclusWithAllCells = true;
% end

thres_reg2 = masterthres;
thres_reg = masterthres;
thres_merge = masterthres;
thres_cap = masterthres;
thres_minsize = 10;
numK1 = 20;

isShowProgress = 1;

% optional override:
if exist('clusParams','var'),
    if ~isempty(clusParams),        
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

%%
autoClusStart = tic;

if isMakeFoxels,
    %% 1. Obtain foxels (functional voxels)
    [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,cIX_reg,isWkmeans,clusParams);%absIX,i_fish);
end

%% 2. Merge foxels to obtain fROI's (functional ROI's)
if isShowProgress,
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

if isShowProgress,
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

%% 
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end

C = FindCentroid_Direct(gIX,M_0(cIX,:));
gIX = HierClus_Direct(C,gIX);
[gIX,numROI] = SqueezeGroupIX(gIX);

%% report timing etc
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