% batch
totaltime = tic; %#ok<NASGU>
isFullData = 1;
data_masterdir = GetCurrentDataDir();

M_stimrange = GetStimRange();

range_fish =  8:18; % range_fish = GetFishRange();

%% custom params here:
% numK1 = 20; 
masterthres = 0.7;

%%
M_regthres = {0.7,0.5};
M_place = {1,2,3,1,4,5,6,7};
M_stimname = {'4x4','PT','OMR','defS','Spt','DF','Lm','Dot'};
%%
for i_count = 1,%%%%%%%%%
    masterthres = M_regthres{i_count};
%     clusParams = struct('merge',masterthres,'cap',masterthres,'reg1',masterthres,...
%         'reg2',masterthres,'minSize',10,'k1',numK1);
    
    for i_stimrange = 1:8,
        if i_stimrange == 1,
            M_stimrange = GetStimRange('5');
        elseif i_stimrange == 2,
            M_stimrange = GetStimRange('P');
        elseif i_stimrange == 3,
            M_stimrange = GetStimRange('O');
        elseif i_stimrange == 4,
            M_stimrange = GetStimRange('M');
        elseif i_stimrange == 5,
            M_stimrange = GetStimRange('S');
        elseif i_stimrange == 6,
            M_stimrange = GetStimRange('D');
        elseif i_stimrange == 7,
            M_stimrange = GetStimRange('L');
        elseif i_stimrange == 8,
            M_stimrange = GetStimRange('Y');
        end

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

            i_ClusGroup = 2;
            i_Cluster = 1;

            % Load cluster data
            [cIX_load,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
                        
            %% partitions for CV
            timelists = getappdata(hfig,'timelists');
            timelists_names = getappdata(hfig,'timelists_names');
            periods = getappdata(hfig,'periods');
            
            M_stim = M_stimrange{i_fish};
            
            timelistsCV_raw = cell(length(M_stim),2);
            timelistsCV = cell(1,2);
            
            for k_stim = 1:length(M_stim), % :3
                i_stim = M_stim(k_stim);
                TL = timelists{i_stim};
                period = periods(i_stim);
                nrep = size(TL,2)/periods(i_stim); % integer
                n = floor(nrep/2);
                if n>0,
                    timelistsCV_raw{k_stim,1} = TL(1:n*period);
                    timelistsCV_raw{k_stim,2} = TL(1+n*period:2*n*period);% before 12/5/16: TL(1+n*period):TL(2*n*period);
                else % for spont, only one period
                    halfperiod = floor(period/2);
                    timelistsCV_raw{k_stim,1} = TL(1:halfperiod);
                    timelistsCV_raw{k_stim,2} = TL(1+halfperiod:2*halfperiod);
                end
            end
            timelistsCV{1} = horzcat(timelistsCV_raw{:,1});
            timelistsCV{2} = horzcat(timelistsCV_raw{:,2});
            assert(length(timelistsCV{1})==length(timelistsCV{2}));
            
            %%
            for k = 1:2,% CV halves
                tIX = timelistsCV{k};
                M = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
                M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');

                % ------custom code here---------
                isWkmeans = true;
                isMakeFoxels = true;
                
                [cIX,gIX] = AutoClustering(cIX_load,gIX,M_0,cIX_load,isWkmeans,[],...
                    isMakeFoxels,masterthres);
                
                % save cluster
                name = ['Auto_',M_stimname{i_stimrange},'_M',num2str(masterthres),'_CV',num2str(k)];
                clusgroupID = 3+k; % 4 and 5
                clusIDoverride = M_place{i_stimrange};
                SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
            end
        end
    end
end
SaveVARwithBackup();
totaltime = toc;
disp(totaltime);

