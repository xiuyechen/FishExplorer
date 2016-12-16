% batch
totaltime = tic; %#ok<NASGU>
isFullData = 1;
data_masterdir = GetCurrentDataDir();

M_stimrange = GetStimRange();

range_fish =  1:18; % range_fish = GetFishRange();

%% custom params here:
% numK1 = 20; 
% masterthres = 0.7;

%%
% M_regthres = {0.7,0.5};
% M_place = {1,2,3,1,4,5,6,7};
% M_stimname = {'4x4','PT','OMR','defS','Spt','DF','Lm','Dot'};

%%

M_GLM_reg = cell(1,length(range_fish));
M_GLM_ctrd = cell(1,length(range_fish));
M_GLM_ctrd_k = cell(1,length(range_fish));
M_GLM_pred_reg = cell(1,length(range_fish));
M_GLM_pred_ctrd = cell(1,length(range_fish));
M_GLM_pred_ctrd_k = cell(1,length(range_fish));

%%
% for i_count = 1,%%%%%%%%%
%     masterthres = M_regthres{i_count};
% %     clusParams = struct('merge',masterthres,'cap',masterthres,'reg1',masterthres,...
% %         'reg2',masterthres,'minSize',10,'k1',numK1);
%     
%     for i_stimrange = 1:8,
%         if i_stimrange == 1,
%             M_stimrange = GetStimRange('5');
%         elseif i_stimrange == 2,
%             M_stimrange = GetStimRange('P');
%         elseif i_stimrange == 3,
%             M_stimrange = GetStimRange('O');
%         elseif i_stimrange == 4,
%             M_stimrange = GetStimRange('M');
%         elseif i_stimrange == 5,
%             M_stimrange = GetStimRange('S');
%         elseif i_stimrange == 6,
%             M_stimrange = GetStimRange('D');
%         elseif i_stimrange == 7,
%             M_stimrange = GetStimRange('L');
%         elseif i_stimrange == 8,
%             M_stimrange = GetStimRange('Y');
%         end

        for i_fishcount = 1:length(range_fish),
            i_fish = range_fish(i_fishcount);
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
            
            %% regressors as feature matrix
%             % CV 1st half: fit model
%             k = 1;
%             tIX = timelistsCV{k};
%             [M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
%             M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
% 
%             fishset = getappdata(hfig,'fishset');
%             [~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
%             [~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
%             Reg = vertcat(regressor_s,regressor_m);
%             regnames = [names_s,names_m]';
%             
%             Data1 = M_0;
%             Reg1 = Reg;
%                 
%             % CV 2nd half: prediction
%             k = 2;
%             tIX = timelistsCV{k};
%             [M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
%             M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
% 
%             fishset = getappdata(hfig,'fishset');
%             [~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
%             [~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
%             Reg = vertcat(regressor_s,regressor_m);
%             regnames = [names_s,names_m]';
%             
%             Data2 = M_0;
%             Reg2 = Reg;
%                 
%              % main code    
% %              [pred,GLM_stat] = GLM_allcells_CV(Data1(1:1000,:),Reg1,Data2(1:1000,:),Reg2);
%             [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1,Reg1,Data2,Reg2);
%             M_GLM_pred_reg{i_fishcount} = GLM_pred_stat;
%             M_GLM_reg{i_fishcount} = GLM_stat;
%             
            %% Autoclus centroid as feature matrix
%             load('C:\Janelia2015\GUI_data\BasicMultiFishInfo.mat');
%             Ctrd = Centroids_Auto07_multi{i_fish}';
            % CV 1st half: fit model
            k = 1;
            tIX = timelistsCV{k};
            % Load cluster data
            i_ClusGroup = 4; % Autocluster based on CV1 --- artifact removed??
            i_Cluster = 1;
            [cIX_load,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
            M = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
            Reg1 = FindCentroid_Direct(gIX,M);
            Data1 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
            
            k = 2;
            tIX = timelistsCV{k};
            M = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
            Reg2 = FindCentroid_Direct(gIX,M);
            Data2 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
            
           % main code    
%            [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1(1:10:end,:),Reg1,Data2(1:10:end,:),Reg2);
% %             [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1,Reg1,Data2,Reg2);%GLM_allcells_CV(Data1(1:1000,:),Reg1,Data2(1:1000,:),Reg2);%
%             M_GLM_pred_ctrd{i_fishcount} = GLM_pred_stat;
%             M_GLM_ctrd{i_fishcount} = GLM_stat; 
            
            % k-means centroid as feature matrix
            numK = size(Reg1,1);
            
            i_ClusGroup = 2; % all cells
            i_Cluster = 1;
            cIX_load = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
            
            k = 1;
            tIX = timelistsCV{k};
            M = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
            [gIX,Ctrd_k] = kmeans(M,numK,'distance','correlation');
            Reg1 = Ctrd_k;
            Data1 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
            
            k = 2;
            tIX = timelistsCV{k};
            M = GetTimeIndexedData_Default_Direct(hfig,cIX_load,tIX);
            Reg2 = FindCentroid_Direct(gIX,M);
            Data2 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');

            % main code
%             [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1(1:10:end,:),Reg1,Data2(1:10:end,:),Reg2);
% %             [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1,Reg1,Data2,Reg2);%GLM_allcells_CV(Data1(1:1000,:),Reg1,Data2(1:1000,:),Reg2);%
%             M_GLM_pred_ctrd_k{i_fishcount} = GLM_pred_stat;
%             M_GLM_ctrd_k{i_fishcount} = GLM_stat; 
            
clusgroupID = 4;
            SaveCluster_Direct(cIX_load,gIX,absIX,i_fish,'kmeans~A0.7',clusgroupID);%,clusIDoverride);
                         
%             % save cluster
%             numK = length(unique(gIX));
%             name = ['k_A0.7_CV1' num2str(numK)];
%             clusgroupID = 3;%+k; % 4 and 5
%             %                 clusIDoverride = M_place{i_stimrange};
%             SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID);%,clusIDoverride);
        end

%     end
% end
% SaveVARwithBackup();
totaltime = toc;
disp(totaltime);

%%
save('GLM_batch_1205_incomplete.mat','M_GLM_pred_reg','M_GLM_reg','M_GLM_pred_ctrd','M_GLM_ctrd','M_GLM_pred_ctrd_k','M_GLM_ctrd_k');

%%
X = M_GLM_pred_ctrd{1};
X(X<-1)=-1;
figure;hist(X,-1:0.05:1)

%%
i_fish = 8;

figure;hold on;
subplot(4,2,1)
X = M_GLM_reg{i_fish}(1:10:end,1);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('reg; fit')

subplot(4,2,2)
X = M_GLM_pred_reg{i_fish}(1:10:end);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('reg; pred')

subplot(423)
X = M_GLM_ctrd{i_fish}(:,1);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('autoclus; fit')

subplot(424)
X = M_GLM_pred_ctrd{i_fish};
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('autoclus; pred')

subplot(425)
X = M_GLM_ctrd_k{i_fish}(:,1);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('kmeans; fit')

subplot(426)
X = M_GLM_pred_ctrd_k{i_fish};
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('kmeans; pred')

%
% subplot(427)
% X = M_GLM_ctrd{i_fish}(:,1)- M_GLM_ctrd_k{i_fish}(:,1);
% X(X<-1)=-1;
% hist(X,-1:0.05:1)
% xlim([-1,1]);ylim([0,5000])
% title('autoclusFit - kmeansFit')
% 
% subplot(428)
% X = M_GLM_pred_ctrd{i_fish}(:,1)- M_GLM_pred_ctrd_k{i_fish}(:,1);
% X(X<-1)=-1;
% hist(X,-1:0.05:1)
% xlim([-1,1]);ylim([0,5000])
% title('autoclusPred - kmeansPred')

subplot(427)
X = M_GLM_ctrd{i_fish}(:,1)- M_GLM_reg{i_fish}(1:10:end,1);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('autoclusFit - regFit')

subplot(428)
X = M_GLM_pred_ctrd{i_fish}(:,1)- M_GLM_pred_reg{i_fish}(1:10:end,1);
X(X<-1)=-1;
hist(X,-1:0.05:1)
xlim([-1,1]);ylim([0,5000])
title('autoclusPred - regPred')
