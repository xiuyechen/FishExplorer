names = fieldnames(VAR);
if strcmp('ClusGroupName',names{1}),
    temp  = struct2cell(VAR);
    VAR = cell2struct(temp, {'ClusFolderName','ClusFolder'});
end