% batch input data setup:
if ~exist('hfig','var'),
    hfig = figure;
end

tic
isFullData = 1;
data_masterdir = GetCurrentDataDir();

M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish = GetFishRange();

%% custom params here:
numK1 = 20;
isWkmeans = true;
isMakeFoxels = false;
% masterthres = 0.5;
% clusParams = struct('merge',masterthres,'cap',masterthres,'reg1',masterthres,...
%     'reg2',masterthres,'minSize',10,'k1',numK1);
M_regthres = {0.7,0.5};

% M_place = {1,2,3,1};
% M_stimname = {'4x4','PT','OMR','defS'};

M_place = {4,5,6,7};
M_stimname = {'Spt','DF','Lm','Dot'};

%%
for i_count = 1:2,
    masterthres = M_regthres{i_count};
    clusParams = struct('merge',masterthres,'cap',masterthres,'reg1',masterthres,...
        'reg2',masterthres,'minSize',10,'k1',numK1);
    
    for i_stim = 1:4,
        if i_stim == 1,
            M_stimrange = GetStimRange('S');
        elseif i_stim == 2,
            M_stimrange = GetStimRange('D');
        elseif i_stim == 3,
            M_stimrange = GetStimRange('L');
        elseif i_stim == 4,
            M_stimrange = GetStimRange('Y');
        end
%         if i_stim == 1,
%             M_stimrange = GetStimRange('5');
%         elseif i_stim == 2,
%             M_stimrange = GetStimRange('P');
%         elseif i_stim == 3,
%             M_stimrange = GetStimRange('O');
%         elseif i_stim == 4,
%             M_stimrange = GetStimRange('M');
%         end
        
        for i = 1:length(range_fish),
            i_fish = range_fish(i);
            disp(i_fish);
            
            % check this loop
            stimrange = M_stimrange{i_fish};
            if isempty(stimrange),
                continue;
            end
            
            % Load fish
            LoadFullFish(hfig,i_fish,isFullData);
            
            %% 1.
            % setup
            absIX = getappdata(hfig,'absIX');
            fishset = M_fishset(i_fish);
            %     i_count = 1;
            i_ClusGroup = 2;
            i_Cluster = 1;

            timelists = getappdata(hfig,'timelists');
            
            % Load cluster data
            [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
            tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
            %     M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
            M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
            
            % ------custom code here---------
            [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
            
            %     % save cluster
%             name = ['Foxels_defStim_reg',num2str(masterthres)];
%             clusgroupID = 5;
%             clusIDoverride = M_place{i_count};
%             SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
            
            % ------custom code here---------
            [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
            
            % save cluster
            name = ['Auto_',M_stimname{i_stim},'_M',num2str(masterthres)];
            clusgroupID = 7+i_count;
            clusIDoverride = M_place{i_stim};
            SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
            
        end
    end
    SaveVARwithBackup();
end
SaveVARwithBackup();
%% 2.
%     % setup
%     i_count = 2;
%     i_ClusGroup = 2;
%     i_Cluster = 2+i_count;
%
%     % Load cluster data
%     [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%     tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%     M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
%
%     % ------custom code here---------
%     [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
%
%     % save cluster
%     name = 'Foxels_half_defStim_reg0.7';
%     clusgroupID = 5;
%     clusIDoverride = i_count;
%     SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%
%     % ------custom code here---------
%     [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
%
%     % save cluster
%     name = 'Auto_half_defStim_Master0.7';
%     clusgroupID = 3;
%     clusIDoverride = i_count;
%     SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);

% end

%% 3.
% for i_fish = 8:18,
%     disp(i_fish);
%
%     % Load fish
%     LoadFullFish(hfig,i_fish,isFullData);
%
%     % setup
%     absIX = getappdata(hfig,'absIX');
%     fishset = M_fishset(i_fish);
%     timelists = getappdata(hfig,'timelists');
%
%     if ismember(i_fish,[8:15,17:18]),
%         % setup
%         i_count = 3;
%         i_ClusGroup = 2;
%         i_Cluster = 2+i_count;
%         stimrange = 1;
%
%         % Load cluster data
%         [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%         tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%         M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
%
%         % ------custom code here---------
%         [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
%
%         % save cluster
%         name = 'Foxels_PT_defStim_reg0.7';
%         clusgroupID = 5;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%
%         % ------custom code here---------
%         [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
%
%         % save cluster
%         name = 'Auto_PT_defStim_Master0.7';
%         clusgroupID = 3;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%
%
%         %%
%         % setup
%         i_count = 4;
%         i_ClusGroup = 2;
%         i_Cluster = 2+i_count;
%         stimrange = 2;
%
%         % Load cluster data
%         [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%         tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%         M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
%
%         % ------custom code here---------
%         [cIX,gIX] = MakeFoxels(cIX,gIX,M_0,isWkmeans,clusParams);
%
%         % save cluster
%         name = 'Foxels_OMR_defStim_reg0.7';
%         clusgroupID = 5;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%
%         % ------custom code here---------
%         [cIX,gIX] = AutoClustering(cIX,gIX,M_0,isWkmeans,clusParams,absIX,i_fish,isMakeFoxels);
%
%         % save cluster
%         name = 'Auto_OMR_defStim_Master0.7';
%         clusgroupID = 3;
%         clusIDoverride = 2+i_count;
%         SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%     end
%
% end
toc

