function [cIX,gIX] = AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans,numK2)
%% set params
% correlation coeff
thres_reg2 = 0.7;
thres_reg = 0.7; 
thres_merge = 0.6; 
thres_cap = 0.6; 
thres_minsize = 10; % number of cells
%%
M = M_0(cIX,:);

%% 1. Obtain 'supervoxels'

%% 1.1. kmeans (2-step)
% step 1:
numK1 = 20;
if isWkmeans,
    disp(['kmeans k = ' num2str(numK1)]);
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

% step 2: divide the above clusters again
% numK2 = 20;
disp(['2nd tier kmeans k = ' num2str(numK2)]);

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

%% 1.2. Regression with the centroid of each cluster
disp('regression with all clusters');
tic
Reg = FindCentroid_Direct(gIX,M);
[cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg,Reg);
gIX = SqueezeGroupIX(gIX);
toc
% clusgroupID = 1;
% SaveCluster_Direct(cIX,gIX,absIX,i_fish,'k20x20_reg',clusgroupID);
%% size thresholding
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
[cIX,gIX] = GrowClustersFromSeedsItr(thres_merge,thres_cap,thres_minsize,thres_reg2,cIX,gIX,M_0);
toc
%%
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
%% Regression again
% C = FindCentroid_Direct(gIX,M_0(cIX,:));
% [cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg2,C);

%% update GUI
C = FindCentroid_Direct(gIX,M_0(cIX,:));
gIX = HierClus_Direct(C,gIX);

clusgroupID = 3;
clusID = SaveCluster_Direct(cIX,gIX,absIX,i_fish,'fromSeed_thres',clusgroupID);
% f.RefreshFigure(hfig);
UpdateClustersGUI_Direct(clusgroupID,clusID,i_fish)

toc
end
