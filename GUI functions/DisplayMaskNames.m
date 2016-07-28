function names_numbered = DisplayMaskNames(Msk_IDs,MASKs)
if Msk_IDs ~= 0,
    names = MASKs.MaskDatabaseNames(Msk_IDs);
    names_numbered = cell(size(names));
    disp('Regions:');
    for i = 1:length(names),
        names_numbered{i} = [num2str(Msk_IDs(i)), ': ', names{i}];
        disp(names_numbered{i});
    end
end
end