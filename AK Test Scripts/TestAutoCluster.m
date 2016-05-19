data_masterdir = GetCurrentDataDir();
range_fish = [2];
%% Load Data from Fish;

%for i = 1:length(range_fish),
    i = 1;
    i_fish = range_fish(i);
    tic
    [F_0,F,N(i),T(i),hfig] = getFullFishData(hfig,i_fish,VAR);
    % it's loading 50% rank right now (see getFullFishData)
    toc
    %end


    %% sample subset of cells
    
    num_select = 1200;%num_cells;
    idx_select = randperm(N(i), num_select);
    
    F_sub = F(idx_select,:);
    
    
    %% Tsne
    tic
    Ftsne =  tsne(F_sub, [], 2);
    toc
    
%%    
    figure(3); hold on;
    plot(Ftsne(:,1),Ftsne(:,2),'.')
    
    
%%  function [cIX,gIX] = AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans,numK2)

%% set params
isWkmeans = true;

% correlation coeff
thres_reg2 = 0.7;
thres_reg = 0.7; 
thres_merge = 0.6; 
thres_cap = 0.6; 
thres_minsize = 4; % number of cells
numK1 = 10;
numK2 = 10;
%%
 M = F_sub;%M_0(cIX,:);
M_0 = F_sub;% for now don't regress to all cells
%% 1. Obtain 'supervoxels'

%% 1.1. kmeans (2-step)
% step 1:

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

%%
figure(3); 
n = round(numK1*1.1);
cmap = hsv(max(1,n));
gscatter(Ftsne(:,1),Ftsne(:,2),gIX,cmap,'.',[],'off');


%%
% step 2: divide the above clusters again

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

%%
figure(4); 
n = round(numK1*1.1);
cmap = hsv(max(1,n));
gscatter(Ftsne(:,1),Ftsne(:,2),gIX,cmap,'.',[],'off');


%% 1.2. Regression with the centroid of each cluster
gIX_old = gIX;
disp('regression with all clusters');
tic
Reg = FindCentroid_Direct(gIX,M);
[cIX,gIX,~] = AllCentroidRegression_direct(M_0,thres_reg,Reg);
gIX = SqueezeGroupIX(gIX);
toc
% clusgroupID = 1;
% SaveCluster_Direct(cIX,gIX,absIX,i_fish,'k20x20_reg',clusgroupID);
%%
figure(5); 
n = round(numK1*numK2*1.1);
cmap = hsv(max(1,n));
gscatter(Ftsne(cIX,1),Ftsne(cIX,2),gIX,cmap,'.',[],'off');

%% size thresholding
U = unique(gIX);
numU = length(U);
gIX_old = gIX;
for i=1:numU,
    if length(find(gIX==U(i)))<thres_minsize/2,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
%[gIX,numU] = SqueezeGroupIX(gIX);
disp(numU);


%% figure(6); 
n = round(numU*1.1);
cmap = hsv(max(1,n));
gscatter(Ftsne(cIX,1),Ftsne(cIX,2),gIX,cmap,'.',[],'off');



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
