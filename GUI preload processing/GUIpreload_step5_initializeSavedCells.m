% initialize VAR (once)
clear all;close all;clc
% clearvars -except 'CellResp' 'const'; clc

data_dir = GetCurrentDataDir();
save_dir = GetCurrentDataDir();
i_fish = 1;

%% load data
disp(['load fish ' num2str(i_fish)]);

filename = ['Data_F' num2str(i_fish) '.mat'];
load(fullfile(data_dir,filename),'data');
% [CellResp,data,dimCR] = LoadFileFromParts(data_dir,filename,'CellResp');

names = fieldnames(data); % cell of strings
for i = 1:length(names),
    eval([names{i} '= data.',names{i} ';']);
end

%% Initialize VAR
VAR = [];
%%
i_ClusGroup = 1;
VAR(i_fish).ClusGroupName{i_ClusGroup} = 'selection';

i_Cluster = 1;
cIX_abs = data.absIX(1:10:end);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '1/10 of validated cells';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;

%
i_ClusGroup = 2;
VAR(i_fish).ClusGroupName{i_ClusGroup} = 'init';

i_Cluster = 1;
cIX_abs = (1:data.numcell_full)';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'all raw cells (ROI''s)';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;

i_Cluster = 2;
cIX_abs = data.absIX;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'validated cells';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;

%%
i_fish = 2;
i_ClusGroup = 1;
i_Cluster = 1;
cIX_abs = (1:data.numcell_full)';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'all raw cells (ROI''s)';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;

i_ClusGroup = 1;
i_Cluster = 2;
cIX_abs = data.absIX;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = 'validated cells';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;

i_ClusGroup = 2;
i_Cluster = 1;
cIX_abs = data.absIX(1:10:end);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).name = '1/10 of validated cells';
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs = cIX_abs;
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).gIX = ones(length(cIX_abs),1);
VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).numK = 1;


%%
save(fullfile(data_dir,'VAR_current.mat'),'VAR');

%%

% varience/std for reps for each cell
% if i_fish==2 || i_fish==3 || i_fish==6,
%     period_real = period/2;
% else
%     period_real = period;
% end
% nrep_real = floor((size(CR,2)-shift)/period_real);
% while period_real*nrep_real+shift>size(CR,2),
%     nrep_real = nrep_real-1;
% end
% CRZ_3D = reshape(CRZ(:,1+shift:period_real*nrep_real+shift),nCells,period_real,[]);
%% updated method, weighing both std between each rep and (summed with) std of 1st half & 2nd half of experiment - 1/8/15
% CRZ = CONST.M_array.CellResp;
% if i_fish==2 || i_fish==3 || i_fish==6,
%     period_real = CONST.M_array.period/2;
% else
%     period_real = CONST.M_array.period;
% end
% CRZ_3D = reshape(CRZ,size(CRZ,1),period_real,[]);
% divide = round(size(CRZ_3D,3)/2);
% CRZ_std1 = std(CRZ_3D(:,:,1:divide),0,3);
% CRZ_std2 = std(CRZ_3D(:,:,divide+1:end),0,3);
% temp1 = mean(CRZ_std1,2);
% temp2 = mean(CRZ_std2,2);
%
% temp12 = horzcat(temp1,temp2);
% temp = mean(temp12,2)+std(temp12,0,2);
% [~,I] = sort(temp);
% M = temp(I);
% figure;plot(M)
%
% figure;imagesc(CRZ(I,:))
%
% nCells = size(CRZ,1);

%% find low variance / stimulus-locked cells
CRZ_std = std(CRZ_3D,0,3);
temp = mean(CRZ_std,2);

% find mean-std thres: 0.5
[~,I] = sort(temp);
M = temp(I);
figure;plot(M)
%%
i_last = length(VAR(i_fish).Class);
M_perc = [0.025,0.1,0.3];
for j = 1:length(M_perc);
    thres = M(round(nCells*M_perc(j)));
    cIX = find(temp<thres);
    i = j+i_last;
    VAR(i_fish).Class(i).round = 0;
    VAR(i_fish).Class(i).name = ['perc < ' num2str(M_perc(j)*100) '%'];
    %     VAR(i_fish).Class(i).notes = ['perc < ' num2str(M_perc(j)*100) '%'];
    VAR(i_fish).Class(i).cIX = cIX;
    VAR(i_fish).Class(i).gIX = ones(length(cIX),1);
    VAR(i_fish).Class(i).numel = length(cIX);
    VAR(i_fish).Class(i).numK = 1;
    VAR(i_fish).Class(i).datatype = 'std';
end

%% shift CR?
% shift = 161;
% nrep=floor(size(CR,2)/period)-1;
%
% skiplist=[];
% IX_rep=setdiff(1:nrep, skiplist);
% IX=zeros(period*length(IX_rep),1);
%
% for i=1:length(IX_rep)
%     IX(period*(i-1)+1:period*i)=period*(IX_rep(i)-1)+1+shift:period*(IX_rep(i))+shift;
% end
%
% % cell_resp=cell_resp(:,1:nrep*period);
% CRA=mean(reshape(CR(:,IX),[nCells,period,length(IX_rep)]),3);
% CRAZ = zscore(CRA')';



