function compareFoldersLR(leftDir,rightDir,newDir,firstNcharMatch)
% This function conveniently assumes that you have matching images in two
% directories ('topDir' and 'bottomDir'), and concatenates the image pairs
% with one on top of the other, and saves it into 'newDir'.
% This allows for easy visual comparison between different sets of plots.
% Option: input a numbber as 'firstNcharMatch' to only match the first N
% characters of the file names.

if ~exist('firstNcharMatch','var')
    isFirstNcharMatch = 0;
else
    isFirstNcharMatch = 1;
end

if ~exist(newDir, 'dir'), mkdir(newDir), end;

% find all file names (that is not a directory, i.e. omit '.' and '..')
listing1 = dir(leftDir);
IX = find(~[listing1.isdir]);

listing2 = dir(rightDir);
IX = find(~[listing2.isdir]);
fileNames2 = cell(1,length(IX));
for i_file = 1:length(IX)
    fileNames2{i_file} = listing2(IX(i_file)).name;
end

for i_file = 1:length(IX)
    % screen files for image files (based on the extension)
    fileName1 = listing1(IX(i_file)).name;
    [~,~,ext] = fileparts(fileName1);
    formatStruct = imformats(ext(2:end)); % first character of 'ext' is '.'
    if ~isempty(formatStruct) % i.e. is matlab-recognized format for imread
        %% check that the corresponding file (with the same name) exists in
        % bottomDir
        
        if isFirstNcharMatch % only match first N characters
            % look in bottomDir for a (first) file that starts with the
            % character string in 'name_check'
            name_check = fileName1(1:firstNcharMatch);
            
            i_file = [];
            M_startIndex = regexp(fileNames2,name_check);
            for i_cell = 1:length(M_startIndex)
                if ~isempty(M_startIndex{i_cell})
                    i_file = i_cell;
                    break;
                end
            end
            if isempty(i_file)
                disp(['Cannot find matching file in bottomDir for file ',name_check,'*']);
                continue;
            end
            fileName2 = fileNames2{i_file};
            
        else
            if exist(fullfile(rightDir,fileName1), 'file')
                fileName2 = fileName1;
            else
                continue;
            end
        end
        
        %% combine images from the 2 directories
        im1 = imread(fullfile(leftDir,fileName1));
        height1 = size(im1,1);
        width1 = size(im1,2);
        im2 = imread(fullfile(rightDir,fileName2));
        height2 = size(im2,1);
        width2 = size(im2,2);
        
        % get new fig size
        h = height1 + height2;
        w = max([width1 width2]);
        figPos = [50 50 w h];
        r = height1/h;
        
        h_limit = 1000;
        if h>h_limit
           r_shrink = h_limit/h;           
           w2 = round(w*r_shrink);
           h2 = h_limit;
           figPos = [50 50 w2 h2];
        end
        w = w2;
        h = h2;
        w_limit = 1916;
        if w>w_limit
            r_shrink = w_limit/w;
            h2 = round(h*r_shrink);
            w2 = w_limit;
            figPos = [50 50 w2 h2];
        end                
        
        % open new figure,  put in subplot:
        f = figure();
        set(f, 'Position', figPos);
        set(f, 'Color', [0 0 0]);
        set(f, 'InvertHardCopy', 'off');
        set(f, 'PaperPositionMode', 'auto');
        subplot('Position',[0 (1-r) 1 r]), imagesc(im1)
        axis('equal'), axis('tight'), axis('off')
        subplot('Position',[0 0 1 (1-r)]), imagesc(im2)
        axis('equal'), axis('tight'), axis('off')
        drawnow();
        
        fn = fullfile(newDir, fileName1);
%         saveas(gcf, fn, 'png');
        print(fn,'-dpng','-r0');
        close(gcf)
    end
end