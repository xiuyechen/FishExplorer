% Scheme incorporating fast dist-test(min(r1,r2)/r12)

%% 1. Pairwise Dist-Test for all clusters from all fish -> candidate clusters
% ~100 clusters each for 18 fish
PairwiseClusterScreen(hfig);
% output: true/false array

%% 2. Threshold pairwise comparison into conserved clusters and save in VAR
ScreenGraphNodes;

