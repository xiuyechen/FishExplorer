% Pool from all fish
clear all;close all;clc

loading_dir = 'C:\Janelia2014';

global VAR;
load('VAR_current.mat','VAR');
%%
numfish = 9;
CONSTs = cell(1,numfish);
for i_fish = 1:numfish,
    disp(['i_fish = ' num2str(i_fish)]);
    
    filename = ['CONST_F' num2str(i_fish) '_fast.mat'];
    [CellResp,const,dimCR] = LoadFileFromParts(loading_dir,filename);
    
    const2 = const;
    const2 = rmfield(const2,{'CellRespAvr','CInfo'});
    CONSTs{i_fish}.const = const2;
%     names = fieldnames(const2); % cell of strings    
    
% additional param
    CONSTs{i_fish}.numcell = length(const.CInfo);
    
%%
    % load selected clusters
    CIX = [];
    for i = 1, % just one for now...
        if length(VAR(i_fish).ClusGroup{1})-i+1>0,
            % copy from VAR
            clus = VAR(i_fish).ClusGroup{1}(i);
            CONSTs{i_fish}.clus(i).cIX = clus.cIX;
            CONSTs{i_fish}.clus(i).gIX = clus.gIX;
            CONSTs{i_fish}.clus(i).name = clus.name;
            CONSTs{i_fish}.clus(i).numK = clus.numK;
            CIX = [CIX; clus.cIX];
        end
    end
    % copy from CONST
    CIX = unique(CIX);
    CONSTs{i_fish}.CIX = CIX;
    CONSTs{i_fish}.CInfo = const.CInfo(CIX); % CIF
    CONSTs{i_fish}.CellRespAvr = const.CellRespAvr(CIX,:);
    CONSTs{i_fish}.CellResp = CellResp(CIX,:);
    
    %% OLD FORMAT without partitioning
%     % load general info
%     load(['CONST_F' num2str(i_fish) '.mat'],'CONST');
%     % CONST: {'ave_stack','anat_yx','anat_yz','CInfo','periods','shift',...
%     %     'CRAZ','CRZt','dshift_fc','FcAvr','Fc','photostate','tlists','datanames'};
% %%
%     % load fields from CONST, with names preserved
%     names = {'anat_yx','anat_yz','anat_zx','periods','shift','dshift','FcAvr','Fc','photostate','tlists','datanames'};
% 
%     for i = 1:length(names),
%         eval(['CONSTs{i_fish}.',names{i},' = CONST.',names{i},';']);
%     end
% 
%     % add one variable
%     CONSTs{i_fish}.numcell = length(CONST.CInfo);
% %%
%     % load selected clusters
%     CONSTs{i_fish}.clus = [];
%     CIX = [];
%     for i = 1, % just one for now...
%         if length(VAR(i_fidsh).ClusGroup{1})-i+1>0,
%             % copy from VAR
%             clus = VAR(i_fish).ClusGroup{1}(i);
%             CONSTs{i_fish}.clus(i).cIX = clus.cIX;
%             CONSTs{i_fish}.clus(i).gIX = clus.gIX;
%             CONSTs{i_fish}.clus(i).name = clus.name;
%             CONSTs{i_fish}.clus(i).numK = clus.numK;
%             CIX = [CIX; clus.cIX];
%         end
%     end
%     % copy from CONST
%     CIX = unique(CIX);
%     CONSTs{i_fish}.CIX = CIX;
%     CONSTs{i_fish}.CInfo = CONST.CInfo(CIX); % CIF
%     CONSTs{i_fish}.CRAZ = CONST.CRAZ(CIX,:);
%     CONSTs{i_fish}.CRZt = CONST.CRZt(CIX,:);
end
%%
temp = whos('CONSTs');
if [temp.bytes]>2*10^9,
    save('CONSTs_current.mat','CONSTs','-v7.3');
else
    save('CONSTs_current.mat','CONSTs');
end

%% optional
timestamp  = datestr(now,'mmddyy_HHMM');
save([loading_dir '\arc mat\CONSTs_' timestamp '.mat'],'CONSTs');
