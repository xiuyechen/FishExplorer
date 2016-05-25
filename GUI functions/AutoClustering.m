function [cIX,gIX] = AutoClusteringAK(cIX,gIX,M_0,isWkmeans,clusParams)
% [cIX,gIX] = AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans,numK2)
% automatically cluster all cells, starting with the currently selected
% cells

%% set params
if exist('clusParams','var'),
    thres_reg2 = clusParams.reg2; % correlation coeff
    thres_reg = clusParams.reg1; % correlation coeff
    thres_merge = clusParams.merge; % correlation coeff
    thres_cap = clusParams.cap; % correlation coeff
    thres_minsize = clusParams.minSize; % number of cells
    numK1 = clusParams.k1;
    numK2 = clusParams.k2;
    
    isShowProgress = 0;
else
    thres_reg2 = 0.5;
    thres_reg = 0.5;
    thres_merge = 0.5;
    thres_cap = 0.5;
    thres_minsize = 10;
    numK1 = 20;
    avrFxsize = 10;
    numK2 = round(length(cIX)/avrFxsize/numK1);
    
    isShowProgress = 1;
end

%%
autoClusStart = tic;
%% 1. Obtain foxels (functional voxels)

%% 1.1. k-means (2-step)
% 1.1.1. first level k-means
k1Start = tic;
if isWkmeans,
    if isShowProgress,
        disp(['kmeans step 1: k = ' num2str(numK1)]);
    end

    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK1 < 10^7 && numK1~=1,
        disp('Replicates = 5');
        gIX = kmeans(M,numK1,'distance','correlation','Replicates',5);
    elseif numel(M)*numK1 < 10^8 && numK1~=1,
        disp('Replicates = 3');
        gIX = kmeans(M,numK1,'distance','correlation','Replicates',3);
    else
        gIX = kmeans(M,numK1,'distance','correlation');
    end
else
    [gIX, numK1] = SqueezeGroupIX(gIX);
end
k1Time = toc(k1Start);

% 1.1.2. divide each of the above clusters again
if isShowProgress,
    disp(['kmeans step 2: k = ' num2str(numK2)]);
end
k2Start = tic;

gIX_old = gIX;
 
% temporarily suppress kmeans warnings
wid = 'stats:kmeans:FailedToConverge';
orig_warn = warning('off',wid);

for i = 1:numK1,
    IX = find(gIX_old == i);
    M_sub = M_0(IX,:);
    rng('default');
    if numK2<length(IX),
        gIX_sub = kmeans(M_sub,numK2,'distance','correlation');
    else
        gIX_sub = kmeans(M_sub,length(IX),'distance','correlation');     
    end
    gIX(IX) = (i-1)*numK2+gIX_sub;
end
% restore warnings
warning(orig_warn);

k2Time = toc(k2Start);

%% 1.2. Regression with the centroid of each cluster
if isShowProgress,
    disp('regression to make foxels');
end
reg1Start = tic;

Reg = FindCentroid_Direct(gIX,M_0);
[cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg,Reg);
gIX = SqueezeGroupIX(gIX);

reg1Time = toc(reg1Start);

%% 1.3 size thresholding
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize/2,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX,numU] = SqueezeGroupIX(gIX);
numFoxels = numU;

%% 2. Merge foxels to obtain fROI's (functional ROI's)
if isShowProgress,
    disp('Merge foxels (iterative)');
end
growStart = tic;
[cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,thres_reg2,cIX,gIX,M_0);
growTime = toc(growStart);

%% 3. Clean-up fROI's
%% 3.1 Regression with the centroid of each cluster (round 2)
if isShowProgress,
    disp('Merge foxels (iterative)');
end
reg2Start = tic;
Reg = FindCentroid_Direct(gIX,M_0(cIX,:));
[cIX,gIX] = AllCentroidRegression_direct(M_0,thres_reg2,Reg);
reg2Time = toc(reg2Start);

%% 3.2. size threshold
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX,numROI] = SqueezeGroupIX(gIX);

%% report timing
autoClusTime = toc(autoClusStart);

disp([' k1:' num2str(k1Time) ...
' k2:' num2str(k2Time) ...
' reg1:' num2str(reg1Time) ... 
' nFox:' num2str(numFoxels) ...
' grow:' num2str(growTime)...
' reg2:' num2str(reg2Time) ...
' nROI:' num2str(numROI) ...
' autoClus:' num2str(autoClusTime) ' sec']);

end