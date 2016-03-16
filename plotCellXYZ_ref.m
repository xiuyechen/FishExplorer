code_dir = 'C:\Users\xiuye\Dropbox\Github\FishExplorer';
addpath(genpath(code_dir));

%%
load('C:\Janelia2015\GUI_data\subject_1\data_full.mat')

%% Import data from text file.
filename = 'C:\Users\xiuye\Downloads\F1_XYZ_scaled_flp_2br_2ref.txt';
% filename = 'C:\Users\xiuye\Downloads\scaled_flipped_transformed\F2_XYZ_scaled_flp_2br_2ref.txt';
delimiter = ' ';
formatSpec = '%f%f%f%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);

VarName1 = dataArray{:, 1};
VarName2 = dataArray{:, 2};
VarName3 = dataArray{:, 3};
FAILED = dataArray{:, 4};

clearvars filename delimiter formatSpec fileID dataArray ans;

CellXYZ_ref = horzcat(VarName1,VarName2,VarName3);

%%
temp = zeros(length(FAILED),1);
for i = 1:length(FAILED),
    temp(i) = length(FAILED{i});
end
    
figure('Position',[100 0 1300 900]);hold on;

subplot(1,2,1)
cIX = find(temp==0);
gIX = round(cIX/1000)+1;
numK = round(cIX(end)/1000)+1;
BasicDrawCellsOnAnatProj(data.CellXYZ,cIX,gIX,numK,data.anat_yx,data.anat_yz);
% BasicDrawCellsOnAnatProj(round(CellXYZ_ref),cIX,gIX,numK,data.anat_yx,data.anat_yz);
title('not FAILED')

subplot(1,2,2)
cIX = find(temp~=0);
gIX = round(cIX/1000)+1;
numK = round(cIX(end)/1000)+1;
BasicDrawCellsOnAnatProj(data.CellXYZ,cIX,gIX,numK,data.anat_yx,data.anat_yz);
% BasicDrawCellsOnAnatProj(round(CellXYZ_ref),cIX,gIX,numK,data.anat_yx,data.anat_yz);
title('FAILED')
