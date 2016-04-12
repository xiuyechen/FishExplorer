% trial using t-test 

range_fish = 8:11;
numFish = length(range_fish);
TF_fishrange = cell(1,numFish);
% also save clusID of matching pairs???
Pairs_AllClusAllFish = []; % each row: {[fish#1, clusID1]}, {[fish#2, clusID2]}

for i_fishnum = 1:numFish,
    % choose reference fish
    i_fish_ref = range_fish(i_fishnum);    
    U_ref = unique(FC{i_fish_ref}.gIX);
    numClus_ref = length(U_ref);
    TF = zeros(numClus_ref,numFish-1);        
    
    % cycle through clusters in reference fish
    for i_ref_clus = 1:numClus_ref,
        refclusID = U_ref(i_ref_clus);
        IX_ref = find(FC{i_fish_ref}.gIX == refclusID);
        Ref_xyz = FC{i_fish_ref}.xyz_norm(IX_ref,:);
        
        % for current cluster, compare with closest cluster in each other fish
        for i_testnum = 1:numFish,
            i_fish_test = range_fish(i_testnum);
            if i_fish_test ~= i_fish_ref,
                U = unique(FC{i_fish_test}.gIX);
                numClus = length(U);
  
                % perform t-test between each pair of clusters
                H = zeros(numClus,3);
                for i = 1:numClus,
                    IX = find(FC{i_fish_test}.gIX==U(i));
                    Test_xyz = FC{i_fish_test}.xyz_norm(IX,:);
                    [h,p] = ttest2(Ref_xyz,Test_xyz,'Vartype','unequal');
                    H(i,:) = h; 
                    % h=0: null hypothesis is true, i.e. clusters match
                end
                
                scores_3D = sum(H,2); % 0: match in all 3 dimensions; 3:in none
                IX = find(scores_3D==0);
                if ~isempty(IX),
                    testclusIX = U(IX);
                    TF(i_ref_clus,i_testnum) = 1;
                    if i_fish_ref<i_fish_test, % the rest is redundant
                        for i = 1:length(testclusIX),
                        Pairs_AllClusAllFish = [Pairs_AllClusAllFish; ...
                            {[i_fish_ref, refclusID]}, {[i_fish_test, testclusIX(i)]}];
                        end
                    end
                else
                    TF(i_ref_clus,i_testnum) = 1-min(scores_3D)/3;
                end
            else % don't compare within same fish
                TF(i_ref_clus,i_testnum) = NaN;
            end
        end        
    end
    TF_fishrange{i_fishnum} = TF;
end

%% summarize conserved clusters as graph
P1 = vertcat(Pairs_AllClusAllFish{:,1});
P2 = vertcat(Pairs_AllClusAllFish{:,2});
Nodes = unique([P1;P2],'rows');
numNodes = length(Nodes);

[~,IX1] = ismember(P1,Nodes,'rows');
[~,IX2] = ismember(P2,Nodes,'rows');

G = graph(IX1,IX2); figure;plot(G)
bins = conncomp(G);

%% Plot anat
anat_yx_norm = getappdata(hfig,'anat_yx_norm');

figure;
hold on;
image(anat_yx_norm)
axis equal
axis ij
axis off

clrmap = hsv(round(max(bins)*1.1));

for i_node = 1:length(Nodes),
   i_fish = Nodes(i_node,1);
   clusIX = Nodes(i_node,2);
   
   GIX = VAR(i_fish).ClusGroup{3}.gIX;
   IX = find(GIX == clusIX);
   CIX_abs = VAR(i_fish).ClusGroup{3}.cIX_abs;
   cIX_abs = CIX_abs(IX);
   LoadFishDataWithoutTS(hfig,i_fish);
   CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
   xyz_norm = CellXYZ_norm(cIX_abs,:);
%    cmap = repmat(clrmap(bins(i_node),:)',1,length(cIX_abs));
   scatter(xyz_norm(:,2),xyz_norm(:,1),2,clrmap(bins(i_node),:));%,20,clrmap(IX,:))
end
