function [cIX,gIX] = AutoClusteringAK(cIX,gIX,M_0,isWkmeans,clusParams)
%% set params
% correlation coeff
thres_reg2 = clusParams.reg2;
thres_reg = clusParams.reg1; 
thres_merge = clusParams.merge; 
thres_cap = clusParams.cap; 
thres_minsize = clusParams.minSize; % number of cells
numK1 = clusParams.k1;
numK2 = clusParams.k2;
%%
%M = M_0(cIX,:);

autoClusStart = tic;
%% 1. Obtain 'supervoxels'

%% 1.1. kmeans (2-step)
% step 1:
k1Start = tic;
if isWkmeans,
    %disp(['kmeans k = ' num2str(numK1)]);
    tic
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
    toc
else
    [gIX, numK1] = SqueezeGroupIX(gIX);
end
k1Time = toc(k1Start);


%% step 2: divide the above clusters again
% numK2 = 20;
k2Start = tic;
%disp(['2nd tier kmeans k = ' num2str(numK2)]);
tic
gIX_old = gIX;
for i = 1:numK1,
    IX = find(gIX_old == i);
    M_sub = M_0(IX,:);
    rng('default');
    if numK2<length(IX),
        [gIX_sub,C] = kmeans(M_sub,numK2,'distance','correlation');
    else
        [gIX_sub,C] = kmeans(M_sub,length(IX),'distance','correlation');
    end
    gIX(IX) = (i-1)*numK2+gIX_sub;
end
k2Time = toc(k2Start);


%% 1.2. Regression with the centroid of each cluster
reg1Start = tic;
%disp('regression with all clusters');

Reg = FindCentroid_Direct(gIX,M_0);
[cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg,Reg);
gIX = SqueezeGroupIX(gIX);

% clusgroupID = 1;
% SaveCluster_Direct(cIX,gIX,absIX,i_fish,'k20x20_reg',clusgroupID);
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

reg1Time = toc(reg1Start);


%% 2.1 Find and Grow Seeds
growStart = tic;
[cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,thres_reg2,cIX,gIX,M_0);
growTime = toc(growStart);

%%
% Regression with the centroid of each cluster, round 2
reg2Start = tic;
Reg = FindCentroid_Direct(gIX,M_0(cIX,:));
[cIX,gIX] = AllCentroidRegression_direct(M_0,thres_reg2,Reg);

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
%disp(numU);

if isempty(gIX),
    errordlg('nothing to display!');
    return;
end

reg2Time = toc(reg2Start);
numClus = numU;


%% update GUI
% C = FindCentroid_Direct(gIX,M_0(cIX,:));
% gIX = HierClus_Direct(C,gIX);
% 
% clusgroupID = 3;
% clusID = SaveCluster_Direct(cIX,gIX,absIX,i_fish,'fromSeed_thres',clusgroupID);
% % f.RefreshFigure(hfig);
% UpdateClustersGUI_Direct(clusgroupID,clusID,i_fish)

% report timing
autoClusTime = toc(autoClusStart);
disp([' k1:' num2str(k1Time) ...
' k2:' num2str(k2Time) ...
' reg1:' num2str(reg1Time) ... 
' nFox:' num2str(numFoxels) ...
' grow:' num2str(growTime)...
' reg2:' num2str(reg2Time) ...
' nClus:' num2str(numClus) ...
' autoClus:' num2str(autoClusTime) ' sec']);


end