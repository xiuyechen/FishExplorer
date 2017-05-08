function [fileNames] = getImageFilesFromFolder(inputDir,namestart_check,fishrange)
%%
listing = dir(inputDir);
IX = find(~[listing.isdir]);

fileNames0 = cell(length(IX),1);
for i_file = 1:length(IX)
    fileNames0{i_file} = listing(IX(i_file)).name;
end

M_num = zeros(1,length(IX));
IX_valid = 1:length(IX);
for i_file = 1:length(IX)
    % screen files for image files (based on the extension)
    fileName = fileNames0{i_file};
    [~,~,ext] = fileparts(fileName);
    formatStruct = imformats(ext(2:end)); % first character of 'ext' is '.'
    if ~isempty(formatStruct) % i.e. is matlab-recognized format for imread
        %% check that the corresponding file (with the same name) exists

        % match first N characters 
        if exist('namestart_check','var')            
            startIndex = regexp(fileName,namestart_check);
            if isempty(startIndex)
                IX_valid(i_file) = 0;
%                 fileNames(i_file) = [];
            end
        end
        
        % only keep files for select fish
        if exist('fishrange','var') 
            num = sscanf(fileName, 'Fish%d*');
            if ~ismember(num,fishrange)
                IX_valid(i_file) = 0;
%                 fileNames(i_file) = [];
            end
            M_num(i_file) = num;
        end
    end
end
IX_valid(IX_valid==0) = [];
[~,B] = sort(M_num(IX_valid));
fileNames = fileNames0(IX_valid(B));
end