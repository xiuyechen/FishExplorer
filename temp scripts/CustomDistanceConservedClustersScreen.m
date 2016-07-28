% AllCentroids{i_fish}.Centroids = C;
%     AllCentroids{i_fish}.XYZn = XYZn;
%     AllCentroids{i_fish}.stim = getappdata(hfig,'stim');
%     AllCentroids{i_fish}.behavior = getappdata(hfig,'behavior');
    
data_dir = GetCurrentDataDir();

load(fullfile(data_dir,'AllCentroids.mat'));

%% pairwise distance
i_fish1 = 8;
i_fish2 = 9;
numClus1 = size(AllCentroids{i_fish1}.Centroids,1);
numClus2 = size(AllCentroids{i_fish2}.Centroids,1);

CostMat = zeros(numClus1,numClus2);
for i = 1:numClus1,%77,%
    i
    Ref_xyz = AllCentroids{i_fish1}.XYZn{i};
    for j = i:numClus2,%32,%
        %%
        Test_xyz = AllCentroids{i_fish2}.XYZn{j};
        
        numCell1 = size(Ref_xyz,1);
        numCell2 = size(Test_xyz,1);
        
        % pairwise cell-distance between clusters
        D1 = pdist(Ref_xyz,'euclidean');
        D2 = pdist(Test_xyz,'euclidean');
        D12 = pdist2(Ref_xyz,Test_xyz,'euclidean');

        [D1_srt,IX] = sort(D1);
        [D2_srt,IX] = sort(D2);
        [D12_srt,IX] = sort(D12(:));

        k = 2;
        D1_avr = mean(D1_srt();
        D2_avr = mean(D2_srt);
        D12_avr = mean(D12_srt);
        
        thres_overlap = 30;
        % first check that cluster is not too spread out
        if min([D1_avr,D2_avr])<100 ... % unit in pixels
                ... % then check for cluster overlap
                && mean(D12_srt(1:round(end/8)))<thres_overlap ...
                    ... % then check that the cluster sizes are not >10x different
                    && length(D1)/length(D2)>1/10 && length(D1)/length(D2)<10,
            score = min([D1_avr,D2_avr])/D12_avr;
        else
            score = 0;
        end
        CostMat(i,j) = score;
        %         if score>1,
        %             % the 2 clusters are similar
        %         end
    end
end

%%
thres_score = prctile(CostMat(:),99); % about 0.5?
IX = find(CostMat(:)>thres_score);
[ind1,ind2] = ind2sub(size(CostMat),IX);
Pairs = zeros(length(IX),4);
Pairs(:,1:3) = [ind1,ind2,CostMat(IX)];

%%
figure;

range_fish = [8,9];
TF = zeros(length(IX),1);
for i = 1:length(IX),
    i
    range_clus = [{Pairs(i,1)},{Pairs(i,2)}]; % good diffuse:[{99},{45}];% sandwiched:[{70},{18}];%
    
    DrawCellsOnAnatProj_MultipleFish(hfig,AllCentroids,range_fish,range_clus,1);
    
    % manually screen clusters: press right mouse button to flag as bad
    waitforbuttonpress
    click = get(gcf,'Selectiontype');
    if strcmp(click,'normal'), % left mouse click
        TF(i) = true;
    elseif strcmp(click,'alt'), % right mouse click
        TF(i) = false;
    end
end

Pairs(:,4) = TF;