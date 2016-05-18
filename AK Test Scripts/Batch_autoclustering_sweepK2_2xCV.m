% batch run full-clustering on all fish

% f.pushbutton_autoclus_Callback

global VAR;

data_masterdir = GetCurrentDataDir();

% range_fish = [5,6,7];
% M_ClusGroup = [2,2,2,2];
% M_Cluster = [1,1,1,1];
range_fish = 2:2;
M_ClusGroup = 1;
M_Cluster = 1;
M_stim = [1];
% M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2];

%%
M_param = 5:5:20%0.3:0.1:0.8;

ParamScores = cell(length(M_param),length(range_fish));
ParamScores_raw = cell(length(M_param),length(range_fish));

for k_param = 1:length(M_param),
    numK2 = M_param(k_param);
    %thres_split = M_param(k_param);
    %setappdata(hfig,'thres_split',thres_split);
    
    for k_fish = 1:length(range_fish),
        i_fish = range_fish(k_fish);
        disp(i_fish);
        LoadFullFish(hfig,i_fish);
        absIX = getappdata(hfig,'absIX');
        
        %% partitions for CV
        timelists = getappdata(hfig,'timelists');
        timelists_names = getappdata(hfig,'timelists_names');
        periods = getappdata(hfig,'periods');
        %if length(periods)>1,
            timelistsCV = cell(length(M_stim),2);
            for k_stim = 1:length(M_stim), % :3
                i_stim = M_stim(k_stim);
                TL = timelists{i_stim};
                period = periods(i_stim);
                nrep = size(TL,2)/periods(i_stim); % integer
                n = floor(nrep/2);
                timelistsCV{k_stim,1} = TL(1):TL(n*period);
                timelistsCV{k_stim,2} = TL(1+n*period):TL(2*n*period);
            end
        %end
        
        %% CV loop: auto-clustering with the partitions
        Score = zeros(length(M_stim),2);
        
        for k_stim = 1:length(M_stim), % :3
            i_stim = M_stim(k_stim);
            NumClus = zeros(1,2);
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
                
                NumClus(k) = length(unique(gIX));
                CIX{k} = cIX;
                GIX{k} = gIX;
            end
            % plot cell-matching figure
            Score(k_stim,1) = HungarianCV(NumClus(1),NumClus(2),CIX{1},CIX{2},GIX{1},GIX{2});
            Score(k_stim,2) = HungarianCV(NumClus(2),NumClus(1),CIX{2},CIX{1},GIX{2},GIX{1});
        end
        
        ParamScores{k_param,k_fish} = mean(mean(Score));
        ParamScores_raw{k_param,k_fish} = Score;
    end
    
end