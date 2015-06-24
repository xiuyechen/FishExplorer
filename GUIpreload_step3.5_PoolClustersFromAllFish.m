% Pool from all fish
clear all;close all;clc

data_dir = 'F:\Janelia2014';

global VAR;
load('VAR_current.mat','VAR');
%%
numfish = 8;
CONSTs = cell(1,numfish);
for i_fish = 1:numfish,
    disp(['i_fish = ' num2str(i_fish)]);
    
    % for new partitioned data, consistant with GUI_050415, but not tested
    fishdir = fullfile(data_dir,['CONST_F' num2str(new_i_fish) '_fast.mat']);
    load(fishdir,'const');
    load(fishdir,'dimCR');
    CellResp = zeros(dimCR);
    num = 0;
    nParts = round(dimCR(1)*dimCR(2)/42000000);
    for i = 1:nParts,
        load(fishdir,['CRZt_' num2str(i)']);
        eval(['len = size(CRZt_' num2str(i) ',1);']);
        eval(['CellResp(num+1:num+len,:) = CRZt_' num2str(i) ';']);
        eval(['CellResp(num+1:num+len,:) = CRZt_' num2str(i) ';']);
        num = num+len;
    end
    names = fieldnames(const); % cell of strings
    for i = 1:length(names),
        % renaming exception!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         if strcmp(names{i},'CRAZ'),
%             setappdata(hfig,'CellRespAvr',const.CRAZ);
%         elseif strcmp(names{i},'photostate'),
%             setappdata(hfig,'stim_full',const.photostate);
        else
            setappdata(hfig,names{i},eval(['const.',names{i}]));
        end
    end
    % add one variable
    CONSTs{i_fish}.numcell = length(const.CInfo);
    
%%
    % load selected clusters
    CONSTs{i_fish}.clus = [];
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
    CONSTs{i_fish}.CRAZ = const.CRAZ(CIX,:);
    CONSTs{i_fish}.CRZt = CellResp(CIX,:);
    
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
%         if length(VAR(i_fish).ClusGroup{1})-i+1>0,
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
clear CONST;
temp = whos('CONSTs');
if [temp.bytes]>2*10^9,
    save('CONSTs_current.mat','CONSTs','-v7.3');
else
    save('CONSTs_current.mat','CONSTs');
end

%% optional
timestamp  = datestr(now,'mmddyy_HHMM');
currentfolder = pwd;
save([currentfolder '\arc mat\CONSTs_' timestamp '.mat'],'CONSTs');
