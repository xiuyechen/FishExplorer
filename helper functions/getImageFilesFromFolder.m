function [fileNames] = getImageFilesFromFolder(inputDir,namestart_check,fishrange)

listing = dir(inputDir);
IX = find(~[listing.isdir]);

fileNames0 = cell(length(IX),1);
for i_file = 1:length(IX)
    fileNames0{i_file} = listing(IX(i_file)).name;
end

fileNames = fileNames0;
for i_file = 1:length(IX)
    % screen files for image files (based on the extension)
    fileName = fileNames0{i_file};
    [~,~,ext] = fileparts(fileName);
    formatStruct = imformats(ext(2:end)); % first character of 'ext' is '.'
    if ~isempty(formatStruct) % i.e. is matlab-recognized format for imread
        %% check that the corresponding file (with the same name) exists in
        % bottomDir
        
        if exist('namestart_check','var') % only match first N characters            
            startIndex = regexp(fileName,namestart_check);
            if isempty(startIndex)
                fileNames(i_file) = [];
            else
                if exist('fishrange','var') % only keep files for select fish
                    str = 
                end
            end
        end
    else
        fileNames(i_file) = [];
    end
end

end