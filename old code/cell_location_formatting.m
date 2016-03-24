dir = 'C:\Users\xiuye\Dropbox\ZBrainOutput\light-sheet data (11 fish)\cell locations';

Collection = cell(1,11);
for i_fish = 1:11,
    i_fish
    fishdir = fullfile(dir,['F',num2str(i_fish),'_cell_info.mat']);
    load(fishdir);
    XY = reshape([cell_info.center],2,[]);
    Z = [cell_info.slice];
    m = vertcat(XY,Z)';
    text = cell(length(m),1);
    for i = 1:length(m),
        text{i} = [num2str(m(i,1)),' ',num2str(m(i,2)),' ',num2str(m(i,3))];
    end
    Collection{i_fish} = text;
end

save(fullfile(dir,'Allfish.mat'),'Collection');

%%

for i_fish = 1:11,
    i_fish
    %%
    fileID = fopen(['F',num2str(i_fish),'_XYZ.txt'],'w');
     
    fishdir = fullfile(dir,['F',num2str(i_fish),'_cell_info.mat']);
    load(fishdir);
    XY = reshape([cell_info.center],2,[]);
    Z = [cell_info.slice];
    m = vertcat(XY,Z)';
    for i = 1:length(m),
       s = [num2str(m(i,2)),' ',num2str(m(i,1)),' ',num2str(m(i,3))];
       fprintf(fileID,[s,'\r\n']);
    end

    fclose(fileID);
    %%
end

%%

for i_fish = 1:11,
    fileID = fopen(['F',num2str(i_fish),'_XYZ.txt'],'w');
    text =Collection{i_fish};
    for i = 1:length(text);
        fprintf(fileID,[text{i},'\r\n']);
    end    
    fclose(fileID);
    %%
end