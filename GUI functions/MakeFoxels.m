function [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,cIX_reg,isWkmeans,clusParams,absIX,i_fish)
% Obtain foxels (functional voxels)

% if isAutoclusWithAllCells,
%     cIX_reg = (1:size(M_0,1))';
% else
%     cIX_reg = cIX;
% end

%% set params
% default:
% thres_reg2 = 0.5;
thres_reg = 0.5;
% thres_merge = 0.5;
% thres_cap = 0.3;
thres_minsize = 10;
numK1 = 20;
avrFxsize = 10;
numK2 = round(length(cIX)/avrFxsize/numK1);

isShowProgress = 1;

% optional override:
if exist('clusParams','var'),
    if ~isempty(clusParams),
        thres_reg = clusParams.reg1; % correlation coeff
        thres_minsize = clusParams.minSize; % number of cells
        numK1 = clusParams.k1;
        
        isShowProgress = 0;
    end
else
    
end

    %% 1. k-means (2-step)
    % 1.1. first level k-means
    k1Start = tic;
    if isWkmeans,
        if isShowProgress,
            disp(['kmeans step 1: k = ' num2str(numK1)]);
        end        
        gIX = Kmeans_Direct(M_0(cIX,:),numK1);
    else
        [gIX, numK1] = SqueezeGroupIX(gIX);
    end
    k1Time = toc(k1Start);
    
    % 1.2. divide each of the above clusters again
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
        M_sub = M_0(cIX(IX),:);
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
    
    %% 2. Regression with the centroid of each cluster
    if isShowProgress,
        disp('regression to make foxels');
    end
    reg0Start = tic;
    
    Reg = FindCentroid_Direct(gIX,M_0(cIX,:));
    [cIX,gIX,numFoxels] = AllCentroidRegression_SizeThres_direct(M_0(cIX_reg,:),cIX_reg,thres_reg,Reg,thres_minsize/2);
    
    reg0Time = toc(reg0Start);
    
%     if exist('absIX','var') && exist('i_fish','var'),
%         clusgroupID = 2;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,'foxels',clusgroupID);
%     end
    
    %% (3. report timing)
    disp(['(Size) nFox:' num2str(numFoxels)]);
    disp(['(Time) k1:' num2str(k1Time) ...
        ' k2:' num2str(k2Time) ...
        ' reg0:' num2str(reg0Time)]);
end