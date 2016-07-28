function [CellXYZ_norm,IX_inval_norm] = GetNormCellCord(i_fish)
if i_fish<=18,
    filename = fullfile('C:\Janelia2015\norm cell coord',['F' num2str(i_fish) '_XYZ_norm.txt']);
    delimiter = ' ';
    formatSpec = '%f%f%f%s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    
    VarName1 = dataArray{:, 1};
    VarName2 = dataArray{:, 2};
    VarName3 = dataArray{:, 3};
    FAILED = dataArray{:, 4};
    
    temp = zeros(length(FAILED),1);
    for i = 1:length(FAILED),
        temp(i) = length(FAILED{i});
    end
    IX_inval_norm = find(temp~=0);
    
    Y = round(VarName1/0.798);
    X = round(VarName2/0.798);
    Z = round(VarName3/2);
    CellXYZ_norm = horzcat(X,Y,Z);
    
else
    CellXYZ_norm = [];
    IX_inval_norm = [];
end
end
