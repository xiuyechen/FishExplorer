% Cross-Validation

%% To rank clusters from 2 different auto-clusterings, for same fish
% maximizing corresponding cell ID's

i_fish = 8;

% first set
i_ClusGroup = 3;
i_Cluster = 1;
Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};

absIX = getappdata(hfig,'absIX');
cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
[~,cIX] = ismember(cIX_abs,absIX);
%     cIX = cIX(find(cIX));
gIX = Cluster(i_Cluster).gIX;

numClus1 = length(unique(gIX));
cIX1 = cIX;
gIX1 = gIX;

% second set
i_ClusGroup = 3;
i_Cluster = 2;
Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};

absIX = getappdata(hfig,'absIX');
cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
[~,cIX] = ismember(cIX_abs,absIX);
%     cIX = cIX(find(cIX));
gIX = Cluster(i_Cluster).gIX;

numClus2 = length(unique(gIX));
cIX2 = cIX;
gIX2 = gIX;

%%
CostMat = zeros(numClus1,numClus2);
for i = 1:numClus1,
    for j = 1:numClus2,
        A = cIX1(gIX1==i);
        B = cIX2(gIX2==j);
        CostMat(i,j) = -length(intersect(A,B));
    end
end

%%
assignment1 = munkres(CostMat);
range = 1:numClus1;
IX = find(assignment1>0);
im1 = -CostMat(range(IX),assignment1(IX));
figure;
imagesc(im1)
colormap(bluewhitered)
axis equal; axis tight
