function [Corr_rows,corr_max] = ConvertCorrToBestRegRows(Corr)
[corr_max,IX_regtype] = max(Corr,[],1);
Corr_rows = zeros(size(Corr));
IX = cell(1,size(Corr,1));
for i_row = 1:size(Corr,1)
    IX{i_row} = find(IX_regtype==i_row);
    Corr_rows(i_row,IX{i_row}) = Corr(i_row,IX{i_row});
end
end