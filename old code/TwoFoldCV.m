function Score = TwoFoldCV(hfig,numK2,M_stim,timelists_names) % need to LoadFullFish first
%% partitions for CV
% timelists = getappdata(hfig,'timelists');
timelists_names = getappdata(hfig,'timelists_names');
% periods = getappdata(hfig,'periods');
% if length(periods)>1,
%     timelistsCV = cell(length(M_stim),2);
%     for k_stim = 1:length(M_stim), % :3
%         i_stim = M_stim(k_stim);
%         TL = timelists{i_stim};
%         period = periods(i_stim);
%         nrep = size(TL,2)/periods(i_stim); % integer
%         n = floor(nrep/2);
%         timelistsCV{k_stim,1} = TL(1):TL(n*period);
%         timelistsCV{k_stim,2} = TL(1+n*period):TL(2*n*period);
%     end
% end

%% CV loop: auto-clustering with the partitions
absIX = getappdata(hfig,'absIX');

Score = zeros(length(M_stim),2);

for k_stim = 1:length(M_stim), % :3
    i_stim = M_stim(k_stim);
    CIX = cell(1,2);
    GIX = cell(1,2);
    for k = 1:2,
        %% Cluster to start auto-clustering
        i_ClusGroup = M_ClusGroup(k_fish);
        i_Cluster = M_Cluster(k_fish);
        ClusGroup = VAR(i_fish).ClusGroup{i_ClusGroup};
        numK = ClusGroup(i_Cluster).numK;
        gIX = ClusGroup(i_Cluster).gIX;
        cIX_abs = ClusGroup(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
        [~,cIX] = ismember(cIX_abs,absIX);
        
        % ~UpdateTimeIndex
        tIX = timelistsCV{k_stim,k};
        M_0 = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX,'isAllCells');
        
        isWkmeans = 1;
        [cIX,gIX] = AutoClustering(cIX,gIX,absIX,i_fish,M_0,isWkmeans,numK2);
        CIX{k} = cIX;
        GIX{k} = gIX;
    end
    % plot cell-matching figure
    Score(k_stim,1) = HungarianCV(CIX{1},CIX{2},GIX{1},GIX{2},timelists_names{i_stim});
    Score(k_stim,2) = HungarianCV(CIX{2},CIX{1},GIX{2},GIX{1},timelists_names{i_stim});
end
end