function [CellResp,data,dimCR] = LoadFileFromParts(data_dir,filename,namestring)
disp('loading file from parts...');
% tic
fishdir = fullfile(data_dir,filename);
load(fishdir,'data');
load(fishdir,'dimCR');
CellResp = zeros(dimCR);
num = 0;
nParts = round(dimCR(1)*dimCR(2)/(10^8));
disp(['in ' num2str(nParts) ' parts']);
for i = 1:nParts,
    disp(num2str(i));
    load(fishdir,[namestring,'_' num2str(i)]);
    eval(['len = size(',namestring,'_' num2str(i) ',1);']);
    eval(['CellResp(num+1:num+len,:) = ',namestring,'_' num2str(i) ';']);
    num = num+len;
end
% toc

end