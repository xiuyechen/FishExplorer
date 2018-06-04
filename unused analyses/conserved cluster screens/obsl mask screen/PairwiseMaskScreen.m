% function PairwiseMaskScreen(hfig,fMaskList_matchZ,fMask_matchZ_ZList);


% FC{i_fish,i_reg} = [];
% FC{i_fish,i_reg}.cIX = cIX;
% FC{i_fish,i_reg}.gIX = gIX;
% FC{i_fish,i_reg}.xyz = xyz;
% FC{i_fish,i_reg}.xyz_avr = xyz_avr;
% FC{i_fish,i_reg}.xyz_norm = xyz_norm;
% FC{i_fish,i_reg}.xyz_norm_avr = xyz_norm_avr;

%% Part 1: from zMasks.

range_fish = 1:18;
numFish = length(range_fish);

% screen ALL pairs within these candidate fMasks
% [instead of "for each zMask, find nodes"]

% for ref fish
%   load fMask_ref range;
%   for each ref mask
%     set fMask_ref;
%     for test fish,   
%         load fMask_test range;    
%         for each test mask,
%             load fMask_test;
%             compute perc overlap for this pair,pool;
% end;end;end;end

% treat ZBrain masks similarly to a ref fish


for i_fishnum = 1:numFish,% cycle through each fish as reference fish

    % choose reference fish
    i_fish_ref = range_fish(i_fishnum); 
    
    % load fMask_ref range
    name = ['fMASKs_',num2str(i_fish),'.mat'];
    load(fullfile(data_masterdir,'fMASKs',name),'fMASKs');
    nfMasks = size(fMASKs.MaskDatabase,2);
    
end


%%
% function PairwiseClusterScreen
range_fish = 2:4;
numFish = length(range_fish);
%% Batch pre-load
numReg = 1;
i_reg = 1;
M_clus = cell(numFish,numReg);
for i_fish = range_fish,
    [cIX,gIX,M_xyz_norm] = GetDefaultClustersFromLoad(hfig,i_fish);
    M_clus{i_fish,i_reg} = [];
%     M_clus{i_fish,i_reg}.cIX = cIX;
    M_clus{i_fish,i_reg}.gIX = gIX;
    M_clus{i_fish,i_reg}.M_xyz_norm = M_xyz_norm;    
end

%% screen ALL pairs (non-directional) between ALL fish
% init
TF_fishrange = cell(1,numFish);
Pairs_AllClusAllFish = []; % each row: {[fish#1, clusID1]}, {[fish#2, clusID2]}

% cycle through each fish as reference fish  
for i_refnum = 1:numFish,
    i_fish_ref = range_fish(i_refnum);
    disp(['i_fish_ref = ' num2str(i_fish_ref)]);
    gIX_ref = M_clus{i_fish_ref,i_reg}.gIX;
    M_xyz_norm_ref = M_clus{i_fish_ref,i_reg}.M_xyz_norm;
    
    % init
    U_ref = unique(gIX_ref);
    numClus_ref = length(U_ref);
    TF_reffish = zeros(numClus_ref,numFish-1);
    
    % cycle through clusters in reference fish
    for i_clus_ref = 1:numClus_ref,
        if mod(i_clus_ref,10)==0,
            disp(['i_clus_ref = ' num2str(i_clus_ref)]);
        end
        clusID_ref = U_ref(i_clus_ref);
        IX = find(gIX_ref == clusID_ref);
        XYZ_ref = M_xyz_norm_ref(IX,:);
        
        % cycle through test fish for given ref fish
        for i_testnum = (i_refnum+1):numFish,% only do ordered pairwise test (other half redundant)
            i_fish_test = range_fish(i_testnum);
            gIX_test = M_clus{i_fish_test,i_reg}.gIX;
            M_xyz_norm_test = M_clus{i_fish_test,i_reg}.M_xyz_norm;
    
            U_test = unique(gIX_test);
            numClus_test = length(U_test);
            scores_test = zeros(numClus_test,1);
            
            % cycle through clusters in test fish
            for i_clus_test = 1:numClus_test,
                clusID_test = U_test(i_clus_test);
                IX = find(gIX_test == clusID_test);
                XYZ_test = M_xyz_norm_test(IX,:);
                
                % compute distance score
                scores_test(i_clus_test,1) = ClusterDistanceTestD12(XYZ_ref,XYZ_test);
            end
            
            % save qualifying pairs for this test fish
            IX_pass = find(scores_test);
            TF_reffish(i_clus_ref,i_testnum) = length(IX_pass);
            if ~isempty(IX_pass),
                clusID_test_pass = U_test(IX_pass);
                for i = 1:length(clusID_test_pass),
                    Pairs_AllClusAllFish = [Pairs_AllClusAllFish; ...
                        {[i_fish_ref, clusID_ref]}, {[i_fish_test, clusID_test_pass(i)]}]; %#ok<AGROW>
                end
            end
        end
    end
    TF_fishrange{i_refnum} = TF_reffish;
end


