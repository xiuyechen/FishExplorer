%% load anat masks (~300)
% original ZBrain anat stack: 1406*621*138
% lightsheet data: e.g. 2048*1188*29
% fMask: cell drawn with radius ~5, smoothen with filter size ~50*50*20
MASKs = getappdata(hfig,'MASKs');

%% 1. make (and store) maps from AutoClus for all fish (~200x18)
MakefMasksForAllFish(hfig);

%% 2. Group1: screen AutoClus with anat masks (~300anat * 200clus * 18fish / 2)
[fMaskList_matchZ,fMask_matchZ_ZList] = ScreenAllfMasksWithzMasks(hfig);

%% 3. Group2: screen AutoClus with regs (~50regs * ~10clus * ~10fish)
%%% (not now, Group3: all AutoClus from all fish, ~200clus * (18*17/2)pairs)
% Batch_poolautoclusters_reg_intoFC;

%% 4. Pool candidate clusters/masks from single fish, and do pairwise comparison
% compare maps within category across fish (~hundreds of candidates only?)
% TtestTrial_loopall
PairwiseMaskScreen(hfig);

%% 5. Pool pairwise comparison into graph to identify top nodes (~100??)
% TtestTrial_loopall

%% 6. save full-size fMask corresponding to best node into annotated collection


%%%%%%%% Scheme incorporating fast dist-test(min(r1,r2)/r12) %%%%%%%%

%% 1. Pairwise Dist-Test for all clusters from all fish -> candidate clusters
% ~100 clusters each for 18 fish
PairwiseClusterScreen(hfig);

%% 2. make (and store) maps for all candidate clusters
MakefMasksForAllFish(hfig);

%% 3. Pairwise Anat-Test for all candidate clusters
PairwiseMaskScreen(hfig);

%% 4. Pool pairwise comparison into graph to identify top nodes (~100??)
ScreenGraphNodes;

%% 5. Screen and label subsets with regs (~50regs * ~10clus * ~10fish)
% Batch_poolautoclusters_reg_intoFC;

%% 6. save full-size fMask corresponding to best node into annotated collection
