% (This is tailored for fishset>1)

clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_dir = GetCurrentDataDir();
% save_dir = GetCurrentDataDir();
i_fish = 11;
%% load data
disp(['load fish ' num2str(i_fish)]);

filename = ['CONST_F' num2str(i_fish) '_fast.mat'];
[CellResp,const,dimCR] = LoadFileFromParts(data_dir,filename);
% fishdir = fullfile(data_dir,filename);
% load(fishdir,'const');

names = fieldnames(const); % cell of strings
for i = 1:length(names),
   eval([names{i} '= const.',names{i} ';']);
end

%% flip
ave_stack = fliplr(ave_stack);
anat_yx = fliplr(anat_yx);
anat_zx = fliplr(anat_zx);

% CInfo = cell info, need to correct indices too
[s1,s2,~] = size(ave_stack);
for i_cell = 1:length(CInfo),
    % fix '.center'
    CInfo(i_cell).center(2) = s2-CInfo(i_cell).center(2)+1;
    % fix '.inds'
    IX = CInfo(i_cell).inds;
    [I,J] = ind2sub([s1,s2],IX);
    J = s2-J+1;
    CInfo(i_cell).inds = sub2ind([s1,s2],I,J);
    % fix '.x_minmax'
    CInfo(i_cell).x_minmax(1) = s2-CInfo(i_cell).x_minmax(1)+1;
    CInfo(i_cell).x_minmax(2) = s2-CInfo(i_cell).x_minmax(2)+1;
end

%%
const = [];
for i = 1:length(names), % use loop to save variables into fields of CONST
    eval(['const.',names{i},'=', names{i},';']);
end
save(fullfile(data_dir,filename),'const','dimCR','-v6');
% custom function:
SaveFileInPartsAppendv6(fullfile(data_dir,filename),CellResp);
