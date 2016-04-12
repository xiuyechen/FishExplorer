% trial using t-test 


range_fish = 8:11;
numFish = length(range_fish);
TF_fishrange = cell(1,numFish);
% also save clusID of matching pairs???
Pairs_AllClusAllFish = []; % each row: fish#1, clusID1, fish#2, clusID2

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
                
                % find testclusID that is closest to current ref-clus in anat-space
                refcoord = mean(Ref.xyz_n,1);
                testcoords = zeros(numClus,3);
                Dist = zeros(numClus,1);
                for i = 1:numClus,
                    IX = find(FC{i_fish_test}.gIX==U(i));
                    testcoords(i,:) = mean(FC{i_fish_test}.xyz_norm(IX,:),1);
                    Dist(i) = pdist([refcoord;testcoords(i,:)]);
                end
                [~,ix] = min(Dist);
                testclusID = U(ix);
                
                % perform t-test between this pair of clusters
                IX = find(FC{i_fish_test}.gIX==testclusID);
                Test_xyz = FC{i_fish_test}.xyz_norm(IX,:);
                [h,p] = ttest2(Ref_xyz,Test_xyz,'Vartype','unequal');
                if max(h) == 0,
                    % clusters are the same, collect
                    TF(i_ref_clus,i_testnum) = 1;
                    if i_fish_ref<i_fish_test,
                        Pairs_AllClusAllFish = [Pairs_AllClusAllFish; ...
                            [i_fish_ref, refclusID, i_fish_test, testclusID]];
                    end
                else
                    TF(i_ref_clus,i_testnum) = 1-sum(h)/3;
                end
            else % don't compare within same fish
                TF(i_ref_clus,i_testnum) = NaN;
            end
        end        
    end
    TF_fishrange{i_fishnum} = TF;
end