clear all; close all; clc
%%
M_dir = GetFishDirectories();
data_masterdir = GetNestedDataDir();
save_masterdir = GetCurrentDataDir();
%%
tiffname = 'C:\Janelia2015\Elavl3-H2BRFP_6dpf_MeanImageOf10Fish.tif';
info = imfinfo(tiffname,'tiff');
nPlanes = length(info);
s1 = info(1).Height;
s2 = info(1).Width;
anat_stack_norm = zeros(s1,s2,nPlanes);
for i=1:nPlanes,
    anat_stack_norm(:,:,i) = imread(tiffname,i);
end
%%
% im = max(anat_stack_norm,[],3);
% out=imNormalize99(im);
% anat_yx_norm = repmat(out,[1 1 3]);
% 
% % y-z view
% im = squeeze(max(anat_stack_norm,[],2));
% out=imNormalize99(im);
% anat_yz_norm = repmat(out,[1 1 3]);
% 
% % x-z view
% im = squeeze(max(anat_stack_norm,[],1));
% out=imNormalize99(im);
% out = flipud(out'); %%%% empirically necessary...
% anat_zx_norm = repmat(out,[1 1 3]);
%%
range_fish = 2:11;

for i_fish = range_fish,
    tic
    disp(['i_fish = ', num2str(i_fish)]);
   
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
%     save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    
    load(fullfile(data_dir,'OptionalInfo.mat'));
%     load(fullfile(save_dir,'data.mat'));
%     absIX = data.absIX;
%     absIX = (1:data.numcell_full)';
    
%     filename = fullfile(save_dir,'TimeSeries.h5');
%     h5create(filename,'/absIX',size(absIX));
%     h5write(filename,'/absIX',absIX);

%     CellXYZ = data.CellXYZ;
%     CellXYZ_ref = data.CellXYZ_ref;
%     
%     load(fullfile(save_dir,'data.mat'));
%     
%     data.CellXYZ = CellXYZ;
%     data.CellXYZ_ref = CellXYZ_ref;
    
%%
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

    %%

    save(fullfile(data_dir,'OptionalInfo.mat'),'Behavior_raw','CellXYZ_norm','IX_inval','numcell_full','IX_inval_norm');%,'anat_stack_norm'

    toc
end

%%
range_fish = 2:11;
for i_fish = range_fish,
    %% load data
    disp(['load fish ' num2str(i_fish) '...']);
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(save_dir,'data_full.mat'),'data'); % struct with many fields
    %%
    data_full = data;
    names = fieldnames(data_full); % cell of strings
    for i = 1:length(names),
        eval([names{i} ' = data_full.' names{i} ';']);
    end
    
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(data_dir,'OptionalInfo.mat'));

    %%
    data = [];
    names = {'periods','timelists_names','stimuluskey_raw','CellXYZ','anat_stack','fpsec',...
        'Behavior_raw','numcell_full','CellXYZ_norm','IX_inval_norm','IX_inval_anat',...
        'anat_yx','anat_yz','anat_zx',...
        'timelists','stim_full','stimAvr','Behavior_full','BehaviorAvr'};
%         'absIX'};
    if length(timelists_names)>1, % M_stimset(i_fish) > 1,
        names = [names,{'stimset'}];
    end
    
    for i = 1:length(names), % use loop to save variables into fields of 'data'
        eval(['data.',names{i},'=', names{i},';']);
    end
    
    save(fullfile(save_dir,'data_full.mat'),'data');
end
%%
filename = fullfile(save_dir,'TimeSeries.h5');
fileattrib(filename,'+w');
h5writeatt(filename,'/','absIX',absIX);
h5disp(filename);

%%
% srcFile = fullfile(matlabroot,'toolbox','matlab','demos','example.h5');
copyfile(filename,'myfile.h5','f');
fileattrib('myfile.h5','+w');
h5writeatt('myfile.h5','/','creation_date',datestr(now));

h5writeatt('myfile.h5','/','absIX',absIX(1:10));

%%
filename = fullfile(save_dir,'TimeSeries_full.h5');
%%
tic
CellResp = h5read(filename,'/CellResp');
CellRespZ = h5read(filename,'/CellRespZ');
CellRespAvr = h5read(filename,'/CellRespAvr');
CellRespAvrZ = h5read(filename,'/CellRespAvrZ');
toc
%%
tic
CellResp_ = CellResp(absIX,:);
CellRespZ_ = CellRespZ(absIX,:);
CellRespAvr_ = CellRespAvr(absIX,:);
CellRespAvrZ_ = CellRespAvrZ(absIX,:);
toc

%%
absIX = (1:data.numcell_full)';

filename = fullfile(save_dir,'TimeSeries_full.h5');
h5create(filename,'/CellResp',size(CellResp),'Datatype','single','ChunkSize',[1000 100]);
h5write(filename,'/CellResp',CellResp);

h5create(filename,'/CellRespZ',size(CellRespZ),'Datatype','single','ChunkSize',[1000 100]);
h5write(filename,'/CellRespZ',CellRespZ);

h5create(filename,'/CellRespAvr',size(CellRespAvr),'Datatype','single','ChunkSize',[1000 100]);
h5write(filename,'/CellRespAvr',CellRespAvr);

h5create(filename,'/CellRespAvrZ',size(CellRespAvrZ),'Datatype','single','ChunkSize',[1000 100]);
h5write(filename,'/CellRespAvrZ',CellRespAvrZ);

h5create(filename,'/absIX',size(absIX));
h5write(filename,'/absIX',absIX);

%% need to delete 'absIX' from data_full.mat

