function [CIX,RegThres] = ChooseTopPercentOfFish(nCells_total,prct_const,Corr_cols)
[Corr_sorted,IX_corr] = sort(Corr_cols,2,'descend');

nCells_target = round(prct_const/100 * nCells_total);

% if 0
%     Msk_IDs = 114; % mask for full hindbrain
%     
%     % isScreenMskFromAllCells
%     cIX = (1:length(absIX))';
%     gIX = ones(size(cIX));
%     [cIX_hb,gIX_hb] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
%     
%     I_hb = ismember(IX_corr,cIX_hb);
%     cum_I_hb1 = cumsum(I_hb(1,:));
%     cum_I_hb2 = cumsum(I_hb(2,:));
%     lastIX{1} = find(cum_I_hb1==nCells_target,1,'first');
%     lastIX{2} = find(cum_I_hb2==nCells_target,1,'first');
% else

%% 
nCol = size(Corr_cols,1);
lastIX = cell(1,nCol);
CIX = cell(1,nCol);
RegThres = cell(1,nCol);
for i_col = 1:nCol
    lastIX{i_col} = nCells_target;
    CIX{i_col} = IX_corr(i_col,1:lastIX{i_col})';
    RegThres{i_col} = Corr_sorted(i_col,lastIX{i_col}); 
end
% end
