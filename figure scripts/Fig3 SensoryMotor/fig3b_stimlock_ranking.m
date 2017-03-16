% fig3B top 1% stimlock
%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load fish
range = 6;%GetFishRange;
for i_fish = range
    ClusterIDs = GetClusterIDs('all');
%     stimrange = 1;
    % [cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
[~,~,M_score] = GetTrialAvrLongTrace(hfig,M);

[score_sort, IX_sort] = sort(M_score);
%% top %
thres_cut = prctile(M_score,2);
thres_min = prctile(M_score,0.1);
ix_end = find(score_sort<thres_cut,1,'last');

IX = IX_sort(1:ix_end);
cIX = cIX_load(IX);
Xrange = [thres_min,thres_cut];
numC = 64;
gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
clrmap = flipud(hot(numC));
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
DrawCellsOnAnat(I);
%% bottom %
thres_cut = prctile(M_score,95);
thres_max = prctile(M_score,99.9);
ix_start = find(score_sort>thres_cut,1,'first');

IX = IX_sort(ix_start:end);
cIX = cIX_load(IX);
Xrange = [thres_cut,thres_max];
numC = 64;
gIX = MapXto1Dcolormap(M_score(IX),Xrange,numC);
clrmap = flipud(hot(numC));
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX,clrmap);
DrawCellsOnAnat(I);
end